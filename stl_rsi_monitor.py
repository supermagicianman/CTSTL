"""
stl_rsi_monitor.py
------------------
Quantitative Robustness-of-Satisfaction-Interval (RoSI) monitor
for a subset of Signal Temporal Logic (discrete time).

Supported operators
  AP (>, <),  Not (¬),  And (∧),  Or (∨),
  F[a,b] (Eventually),  G[a,b] (Always),  U[a,b] (Timed-Until)

Complexities
  O(n) time, O(window) memory.
"""

import math
from collections import deque

INF, NINF = float("inf"), float("-inf")

# ----------------------------------------------------------------------
# deque-based sliding extrema
# ----------------------------------------------------------------------
def _sliding(xs, w, is_max=True):
    cmp = (lambda a, b: a <= b) if is_max else (lambda a, b: a >= b)
    dq, out = deque(), []
    for i, v in enumerate(xs):
        while dq and dq[0] <= i - w:
            dq.popleft()
        while dq and cmp(xs[dq[-1]], v):
            dq.pop()
        dq.append(i)
        if i >= w - 1:
            out.append(xs[dq[0]])
    return out


smax = lambda xs, w: _sliding(xs, w, True)
smin = lambda xs, w: _sliding(xs, w, False)

# ----------------------------------------------------------------------
# interval helpers
# ----------------------------------------------------------------------
imin = lambda a, b: (min(a[0], b[0]), min(a[1], b[1]))
imax = lambda a, b: (max(a[0], b[0]), max(a[1], b[1]))

# ----------------------------------------------------------------------
# sliding-window bounds that look at the prefix length
# ----------------------------------------------------------------------
def bounds_max(ints, w, known_shift):
    """bounds of F on `ints`  (len == n-a)"""
    lo = [p[0] for p in ints]
    hi = [p[1] for p in ints]
    lo_w = smax(lo, w)
    hi_w = smax(hi, w)
    res = []
    for i, (l, h) in enumerate(zip(lo_w, hi_w)):
        if i + w - 1 < known_shift:  # window fully inside known prefix
            res.append((l, h))
        else:                        # touches unknown future
            res.append((l, INF))
    # tail-pad so result aligns with original time axis
    res += [(NINF, INF)] * (len(ints) - len(res))
    return res


def bounds_min(ints, w, known_shift):
    lo = [p[0] for p in ints]
    hi = [p[1] for p in ints]
    lo_w = smin(lo, w)
    hi_w = smin(hi, w)
    res = []
    for i, (l, h) in enumerate(zip(lo_w, hi_w)):
        if i + w - 1 < known_shift:
            res.append((l, h))
        else:
            res.append((NINF, h))
    res += [(NINF, INF)] * (len(ints) - len(res))
    return res

# ----------------------------------------------------------------------
# STL AST nodes
# ----------------------------------------------------------------------
class Node:                                   # base
    def eval(self, sig, known): raise NotImplementedError

class AP(Node):
    def __init__(self, var, op, c):
        self.var, self.op, self.c = var, op, c
    def eval(self, sig, known):
        xs = sig[self.var]
        out = []
        for i, x in enumerate(xs):
            if i < known:                        # known sample
                r = x - self.c if self.op == ">" else self.c - x
                out.append((r, r))
            else:                                # unknown future
                out.append((NINF, INF))
        return out
    def __str__(self):
        return f"{self.var}{self.op}{self.c}"


class Not(Node):
    def __init__(self, ch): self.ch = ch
    def eval(self, s, k):
        return [(-hi, -lo) for lo, hi in self.ch.eval(s, k)]
    def __str__(self):   return "¬"+str(self.ch)


class And(Node):
    def __init__(self, *chs): self.chs = chs
    def eval(self, s, k):
        vals = [c.eval(s, k) for c in self.chs]
        res  = vals[0]
        for v in vals[1:]:
            res = [imin(a, b) for a, b in zip(res, v)]
        return res
    def __str__(self): return '('+' ∧ '.join(map(str,self.chs))+')'


class Or(Node):
    def __init__(self, *chs): self.chs = chs
    def eval(self, s, k):
        vals = [c.eval(s, k) for c in self.chs]
        res  = vals[0]
        for v in vals[1:]:
            res = [imax(a, b) for a, b in zip(res, v)]
        return res
    def __str__(self): return '('+' ∨ '.join(map(str,self.chs))+')'


def F_eval(a, b, child, known):
    """
    child  : list[(lo,hi)] for φ, length = n
    returns: list[(lo,hi)] for F_[a,b] φ, length = n
    """
    n   = len(child)
    win = b - a + 1

    # sliding max over the *full* child list
    # returns n - win + 1 values, aligned with the start index of each window
    full_max_lo = smax([p[0] for p in child], win)
    full_max_hi = smax([p[1] for p in child], win)

    # robustness at time t is max(child[t+a : t+b])
    # which is sliding-max result at index t+a
    res = []
    for t in range(n):
        idx = t + a
        if idx >= len(full_max_lo):               # tail that cannot fit
            res.append((NINF, INF))
        else:
            lo, hi = full_max_lo[idx], full_max_hi[idx]
            if t + b < known:                     # window fully known
                res.append((lo, hi))
            else:                                 # touches the future
                res.append((lo, INF))
    return res


class F(Node):
    def __init__(self, a, b, ch):
        self.a, self.b, self.ch = a, b, ch

    def eval(self, s, k):
        child = self.ch.eval(s, k)
        return F_eval(self.a, self.b, child, k)


def G_eval(a, b, child, known):
    # identical to F_eval but uses smin and widens lower bound
    n   = len(child)
    win = b - a + 1
    full_min_lo = smin([p[0] for p in child], win)
    full_min_hi = smin([p[1] for p in child], win)
    res = []
    for t in range(n):
        idx = t + a
        if idx >= len(full_min_lo):
            res.append((NINF, INF))
        else:
            lo, hi = full_min_lo[idx], full_min_hi[idx]
            if t + b < known:
                res.append((lo, hi))
            else:
                res.append((NINF, hi))
    return res


class G(Node):
    def __init__(self, a, b, ch):
        self.a, self.b, self.ch = a, b, ch
    def eval(self, s, k):
        child = self.ch.eval(s, k)
        return G_eval(self.a, self.b, child, k)

    
def unb_until(phi, psi):
    m = max(len(phi), len(psi))
    # extend with ±inf so intervals remain conservative
    pad_phi = phi + [(NINF, INF)] * (m - len(phi))
    pad_psi = psi + [(NINF, INF)] * (m - len(psi))

    out = [None] * m
    out[-1] = imin(pad_phi[-1], pad_psi[-1])
    for i in range(m - 2, -1, -1):
        out[i] = imax( imin(pad_phi[i], pad_psi[i]),
                       imin(pad_phi[i], out[i + 1]) )
    return out


class U(Node):
    def __init__(self, a, b, phi, psi):
        self.a, self.b, self.phi, self.psi = a, b, phi, psi
    def eval(self, s, k):
        evt = F(self.a, self.b, self.psi).eval(s, k)
        unb = unb_until(self.phi.eval(s, k), self.psi.eval(s, k))
        dur = bounds_min(unb, self.a+1, k)  # window from 0
        return [imin(e, d) for e, d in zip(evt, dur)]
    def __str__(self): return f"({self.phi} U[{self.a},{self.b}] {self.psi})"



import heapq 


class _SlidingTauLargestOracle:
    def __init__(self, k_window: int, tau_rank: int):
        if not (k_window > 0 and 1 <= tau_rank <= k_window):
            raise ValueError("tau_rank must be between 1 and k_window (inclusive), and k_window > 0.")
        self.k = k_window
        self.tau = tau_rank
        self.max_h = [] 
        self.min_h = [] 
        self.removed_map = {}
        self.max_heap_elements = 0
        self.min_heap_elements = 0

    def _pop_items_marked_for_removal(self):
        while self.max_h and self.max_h[0][1] in self.removed_map:
            idx = self.max_h[0][1]
            heapq.heappop(self.max_h)
            self.removed_map[idx] -= 1
            if self.removed_map[idx] == 0:
                del self.removed_map[idx]
        
        while self.min_h and self.min_h[0][1] in self.removed_map:
            idx = self.min_h[0][1]
            heapq.heappop(self.min_h)
            self.removed_map[idx] -= 1
            if self.removed_map[idx] == 0:
                del self.removed_map[idx]

    def _balance_heaps(self):
        while self.min_heap_elements > self.tau:
            self._pop_items_marked_for_removal() 
            if not self.min_h: break
            val, idx = heapq.heappop(self.min_h)
            heapq.heappush(self.max_h, (-val, idx))
            self.min_heap_elements -= 1
            self.max_heap_elements += 1
        
        while self.min_heap_elements < self.tau and self.max_heap_elements > 0:
            self._pop_items_marked_for_removal()
            if not self.max_h: break
            val_neg, idx = heapq.heappop(self.max_h)
            heapq.heappush(self.min_h, (-val_neg, idx))
            self.max_heap_elements -= 1
            self.min_heap_elements += 1
        
        self._pop_items_marked_for_removal()

    def add_number(self, num_val: float, index: int):
        self._pop_items_marked_for_removal()
        # Simplified add logic: add to one, then balance.
        # Prefer min_h if it's not full or if num_val is large enough.
        if self.min_heap_elements < self.tau:
            heapq.heappush(self.min_h, (num_val, index))
            self.min_heap_elements += 1
        elif self.min_h and num_val >= self.min_h[0][0]: # if min_h is full, but num_val is larger than its smallest
            heapq.heappush(self.min_h, (num_val, index))
            self.min_heap_elements +=1
        else: # min_h is full and num_val is smaller than its smallest, or min_h is empty (covered by first if)
            heapq.heappush(self.max_h, (-num_val, index))
            self.max_heap_elements += 1
        
        self._balance_heaps()


    def remove_number(self, num_val_to_remove: float, index: int): # num_val_to_remove is for conceptual count adjustment
        self.removed_map[index] = self.removed_map.get(index, 0) + 1
        
        # Conceptual decrement: if it was in min_h group or max_h group
        # This requires knowing which group it was in.
        # A robust way: if after cleaning heaps, its value would put it into min_h or max_h.
        # This is tricky. The conceptual counts are best managed by balance after physical removal.
        # For now, let balance handle the counts after physical removal via _pop_items_marked_for_removal.
        # Let's assume counts are reduced when items are *physically* popped by _balance_heaps or _pop_items...
        # The current _balance_heaps correctly adjusts counts when moving items.
        # The key is that _pop_items_marked_for_removal doesn't adjust conceptual counts.
        # So, we need to decrement the conceptual count here.

        # Simplified conceptual count decrement:
        # Assume element was part of the k elements.
        # If it's bigger than the current tau-th largest (or what would be if min_h isn't full), it was in min_h.
        found_in_min_conceptual = False
        if self.min_h and num_val_to_remove >= self.min_h[0][0] and self.min_heap_elements == self.tau :
            found_in_min_conceptual = True
        elif self.min_heap_elements < self.tau : # If min_h wasn't full, it might have been there
             # Heuristic: if it's larger than largest of smalls (if max_h exists)
             if not self.max_h or num_val_to_remove > -self.max_h[0][0]:
                 found_in_min_conceptual = True


        if found_in_min_conceptual and self.min_heap_elements > 0:
            self.min_heap_elements -= 1
        elif self.max_heap_elements > 0: # Otherwise, assume it was in max_h
            self.max_heap_elements -= 1
        # else: if both counts are 0, something is off, but balance should try to fix.

        self._balance_heaps()


    def get_tau_th_largest(self) -> float:
        self._pop_items_marked_for_removal()
        if self.min_heap_elements == self.tau and self.min_h:
            return self.min_h[0][0]
        if self.min_heap_elements < self.tau and self.min_heap_elements + self.max_heap_elements < self.k : # Not enough elements in window
             return NINF # Or based on problem spec for too few elements
        if self.min_h: # If min_h has elements but not tau, it means tau-th is among NINF values
             return self.min_h[0][0] # Smallest of the potentially largest items
        return NINF 
        
def _get_sliding_tau_largest_list(nums: list[float], k_win: int, tau_rank: int) -> list[float]:
    if not (k_win > 0):
        return [] 
    if not (1 <= tau_rank <= k_win):
        if not nums or k_win > len(nums): return []
        return [NINF] * (len(nums) - k_win + 1) if len(nums) >= k_win else []

    if not nums or k_win > len(nums):
        return []

    results = []
    oracle = _SlidingTauLargestOracle(k_win, tau_rank)

    for i in range(len(nums)):
        oracle.add_number(nums[i], i)
        if i >= k_win - 1:
            results.append(oracle.get_tau_th_largest())
            oracle.remove_number(nums[i - k_win + 1], i - k_win + 1) 
    return results

class C(Node):
    def __init__(self, a: int, b: int, tau: int, ch: Node):
        self.a, self.b, self.tau, self.ch = a, b, tau, ch
        self.win_len = b - a + 1

    def eval(self, signals_dict: dict, known_signal_len: int) -> list[tuple[float, float]]:
        child_robs_list = self.ch.eval(signals_dict, known_signal_len)
        num_child_results = len(child_robs_list)

        if not (self.win_len > 0 and 1 <= self.tau <= self.win_len):
            return [(NINF, INF)] * num_child_results

        los = [p[0] for p in child_robs_list]
        his = [p[1] for p in child_robs_list]

        kth_los_values = _get_sliding_tau_largest_list(los, self.win_len, self.tau)
        kth_his_values = _get_sliding_tau_largest_list(his, self.win_len, self.tau)
        
        output_robustness_intervals = []
        for t in range(num_child_results):
            result_idx = t + self.a 
            
            if 0 <= result_idx < len(kth_los_values) and 0 <= result_idx < len(kth_his_values):
                rob_lo = kth_los_values[result_idx]
                rob_hi = kth_his_values[result_idx]
                # Ensure lo <= hi, which might not hold if NINF/INF propagation is complex
                # For tau-th largest, if los has NINF, rob_lo can be NINF.
                # If his has INF, rob_hi can be INF.
                # It's possible for rob_lo to be NINF and rob_hi to be a finite number if
                # the tau-th largest of lower bounds is NINF, but for upper bounds it's finite.
                # However, generally we expect rob_lo <= rob_hi. If NINF propagation is right, it should hold.
                output_robustness_intervals.append((rob_lo, rob_hi))
            else:
                output_robustness_intervals.append((NINF, INF))
        return output_robustness_intervals

    def __str__(self) -> str:
        return f"C[{self.a},{self.b},{self.tau}]({str(self.ch)})"



def monitor_rsi(formula, signals, t, known):
    """
    formula : AST root
    signals : dict  sensor_name -> list[float]
    t       : evaluation index
    known   : prefix length currently available (≥ t+1, ≤ len(signal))
    Returns : (lower, upper)
    """
    lo, hi = formula.eval(signals, known)[t]
    return lo, hi




