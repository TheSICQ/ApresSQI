# AprèsSQI-POC

Implementation of the SQIsign signature algorithm, using field extensions when signing, allowing for lightning fast verification. 

### Example usage
```python
from ApresSQI import SQIsign

# Signer
msg = "Kani came in like a wrecking ball"
signer = SQIsign()
pk = signer.KeyGen()
signer.sk.export('priv.key')

sigma = signer.Sign(msg)

# Verifier
verifier = SQIsign()
verified = verifier.verify(msg, sigma, pk)
print('Signature was ' + ('CORRECT' if verified else 'WRONG'))
```

Or, if you just want to test a ready-made example, you can run
```bash
sage example.py
```
which allows for signing with the four parameters 'toy', 'NIST', '7-block', '4-block', as well as signing both with and without seeds, or uncompressed.

(tested with SageMath version 9.8)

### Acknowledgements
The code is based on the [documentation](http://sqisign.org/spec/sqisign-20230601.pdf) from the NIST submission of SQIsign. Further, it uses many ideas introduced in [Deuring for the People](https://eprint.iacr.org/2023/106), and some code is also taken from the [code](https://github.com/friends-of-quaternions/deuring), made by Jonathan Komada Eriksen, Lorenz Panny, Jana Sotáková and Mattia Veroni. Finally, some [code](https://github.com/LearningToSQI/SQISign-SageMath) is also taken from the project [Learning to SQI](https://learningtosqi.github.io/), made by Maria Corte-Real Santos, Jonathan Komada Eriksen, Michael Meyer and Giacomo Pope.