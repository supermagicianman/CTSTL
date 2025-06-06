# CT-STL: Cumulative-Time Signal Temporal Logic Monitor

We provide a online monitor for CT-STL formula and discrete-time signal with unit time step 1.

The quantitative semantics of the cumulative operator C^tau_[a,b] is the tau-th largest value on time interval [a,b]. Other operators are the same as in STL.

The monitor code are defined in stl_rsi_monitor.py. The online monitoring returns Robustness Satisfaction Interval (RoSI), which is a pair including lower and upper bound of the robustness over a partial signal.

The AP (Artificial Pancreas) and microgrid are the case studies in our paper: Cumulative-Time Signal Temporal Logic by Hongkai, et al.

## Usage

The monitor example is in monitor_tutorial.ipynb.

```python
from stl_rsi_monitor import *

# define input signal e.g., {"x": [5,10,12,-3,4,-2], "y":[-3,9,4,10,-3,5], ...}
signals = {'x': [2, -1, 7, 10, -5, 15, 8, -2]}
print(signals)

# define CTSTL formula
# Atomic proposition: AP(signal_name:str, op:str, c:int))     x>c; x<c
# Negation: Not(phi:Node)                      Not phi
# And: And(phi1:Node, phi2:Node, ...)           phi1 And phi2 And ...
# Or: Or(phi1:Node, phi2:Node, ...)           phi1 Or phi2 Or ...
# Always: G(a:int, b:int, phi:Node)           G_[a,b] phi                  
# Eventually: F(a:int, b:int, phi:Node)           F_[a,b] phi
# Until: U(a:int, b:int, phi1:Node, phi2:Node)           phi1 U_[a,b] phi2

# G_[0,2] C^3_[1,5] x>0
formula = G(0,2,
           C(1,5,3,
            AP("x",">",0)))


# calculate rosi at t=0, t should be less than signal length
t=0  

# known is the length of the signal prefix which are already known, known should be less than the signal length
known = len(signals['x']) 

# rosi = monitor_rsi(formula:Node, signals, t, known)
lo, hi = monitor_rsi(formula, signals, t=t, known=known)
print(f"t={t:2d}, prefix={known:2d}  →  RoSI = [{lo: .1f}, {hi: .1f}]")

# assume partial signal is known
lo, hi = monitor_rsi(formula, signals, t=t, known=known-1)
print(f"t={t:2d}, prefix={known-1:2d}  →  RoSI = [{lo: .1f}, {hi: .1f}]")

lo, hi = monitor_rsi(formula, signals, t=t, known=known-2)
print(f"t={t:2d}, prefix={known-2:2d}  →  RoSI = [{lo: .1f}, {hi: .1f}]")

lo, hi = monitor_rsi(formula, signals, t=t, known=known-3)
print(f"t={t:2d}, prefix={known-3:2d}  →  RoSI = [{lo: .1f}, {hi: .1f}]")

# output:
# {'x': [2, -1, 7, 10, -5, 15, 8, -2]}
# t= 0, prefix= 8  →  RoSI = [ 7.0,  7.0]
# t= 0, prefix= 7  →  RoSI = [ 7.0,  7.0]
# t= 0, prefix= 6  →  RoSI = [-5.0,  7.0]
# t= 0, prefix= 5  →  RoSI = [-inf,  7.0]
```
