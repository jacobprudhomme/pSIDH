# pSIDH
Implementation of pSIDH post-quantum cryptographic key exchange described in https://eprint.iacr.org/2021/1600.

## Prerequisites

## Getting Started

## Using the Library

## Todo

[ ] Wait for Dr. Leroux to give the precise meaning of Z + DO

[ ] Wait for Dr. Leroux's answer on the nature of the alpha variable found in KeyExchange

[ ] Implement multidimensional DLP for use in SuborderEvaluation

[ ] Implement efficient way of obtaining the integer T we need in KeyExchange

[ ] Verify that RandomEquivalentPrimeIdeal finds an ideal of relatively close norm to the original. If so, implement a new way of generating an ideal of norm D. If not, go back to the drawing board

[ ] Include all assertions found in the pseudocode

[ ] Wrap key exchange in class to hold global parameters we are currently generating externally. Use constructor to set these

[ ] Allow public-key validation to be omitted

[ ] Verify that suborder representations are correctly represented in code
