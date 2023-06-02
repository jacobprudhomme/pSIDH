# pSIDH
Implementation of pSIDH post-quantum cryptographic key exchange described in https://eprint.iacr.org/2021/1600, for my semester research project in EPFL's [LASEC](https://lasec.epfl.ch) lab.

## Prerequisites

Sage release 9.8 is required to run this code. Compatibility is not guaranteed with newer versions.

## Getting Started

Clone the repo with the `git clone --recurse-submodules` command. If the repo has already been cloned, enter its directory and run `git submodule update --init`. This fetches the other Sage libraries we rely on.

## Using the Library

Once the library is in a fully working state, the process to use it should look something like this: import the pSIDH class and initialize it with one of the possible parameters p. Then, we can run the key generation and key exchange using this object:

```python
from pSIDH import pSIDH
import pSIDH.params.p_toy as params  # Or p6983

inst = pSIDH(p)

alice_public_key, alice_secret_key, D = inst.KeyGeneration()
bob_public_key, bob_secret_key, D_prime = inst.KeyGeneration()

alice_shared_secret = inst.KeyExchange(alice_secret_key, D_prime, ...bob_public_key) # By default, does public-key validation
bob_shared_secret = inst.KeyExchange(bob_secret_key, D, ...alice_public_key, validate=False) # Do not validate public-key

assert alice_shared_secret == bob_shared_secret, 'Key exchange failed: the shared secrets are not the same'
```

## Todo

- [ ] Wait for Dr. Leroux to give the precise meaning of Z + DO

- [ ] Wait for Dr. Leroux's answer on the nature of the alpha variable found in KeyExchange

- [ ] Implement multidimensional DLP for use in SuborderEvaluation

- [ ] Implement efficient way of obtaining the integer T we need in KeyExchange

- [ ] Verify that RandomEquivalentPrimeIdeal finds an ideal of relatively close norm to the original. If so, implement a new way of generating an ideal of norm D. If not, go back to the drawing board

- [ ] Include all assertions found in the pseudocode

- [ ] Wrap key exchange in class to hold global parameters we are currently generating externally. Use constructor to set these

- [ ] Allow public-key validation to be omitted

- [ ] Verify that suborder representations are correctly represented in code
