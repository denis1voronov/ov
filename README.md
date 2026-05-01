# Complexity estimations
We use cryptographic-estimators package (https://github.com/Crypto-TII/CryptographicEstimators), licensed under Apache-2.0. To run attack estimator, execute the following:

```bash
git clone https://github.com/denis1voronov/ov
cd ov
python -m venv .venv

# Windows
.venv\Scripts\activate

# Linux/macOS
source .venv/bin/activate

pip install -r requirements.txt
python complexity-estimations/estimations.py
```

# Algorithm 1 implementation
To run attack experiments, make sure to have SageMath 10.7 available. Additionally, msolve 0.9.4 library (https://msolve.lip6.fr/) is required. It is used for efficient computation of quadratic systems solutions over prime fields. Over extension fields, native Sage is used. Then, run

```bash
sage implementation/attack.sage
```
