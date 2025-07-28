# SGSA: Stage GCD Sieving Algorithm
**Version**: 1.00  
**Language**: PARI/GP

This script implements the Stage GCD Sieving Algorithm (SGSA), a method for efficiently sieving prime candidates by computing greatest common divisors (GCDs) with precomputed prime products for each stage up to a predetermined limit. It filters out composite numbers before performing more computationally expensive primality tests, such as the Fermat primality test. 

## Overview

The algorithm can be used in two modes: sequential mode and batch mode.

Sequential mode is suitable for applying SGSA to a single number or to a group of large numbers, processed one at a time. 

Batch mode is designed to process a single batch containing many numbers simultaneously, provided that sufficient hardware resources are available. It outperforms sequential mode in most scenarios.

SGSA can be applied to any number without restrictions. While it may not be as efficient as specialized sieving methods designed for numbers generated in a specific manner, for general number inputs, our tests show that it offers significant performance advantages over methods such as trial factoring, Pollard’s rho, and ECM in certain scenarios. 

Additionally, the script includes two functions for generating prime candidates using an adaptation of Pocklington’s theorem.

Further details on SGSA, its comparison with other methods, and the candidate generation process based on Pocklington’s theorem can be found in the following paper: https://doi.org/10.36227/techrxiv.175373801.16133264/v1

This script may be particularly useful for researchers in number theory and anyone interested in optimizing prime search algorithms.

## Requirements
### PARI/GP 
PARI/GP 2.15.5 (x86-64/GMP-6.1.2 kernel) is required. Other versions might work, but compatibility is not guaranteed.

Website: https://pari.math.u-bordeaux.fr/

### RAM
Sufficient memory is required for handling large prime product computations.

The required PARI/GP available RAM depends on the size of the sieve stages: 

* For products of 10⁷ primes, `make_products(10,7,max_stage)`, around 2.8GB of available RAM in PARI/GP is required. To set this, use: `default(parisize, 3000000000)` 
* For products of 10⁸ primes, `make_products(10,8,max_stage)`, around 28GB of available RAM in PARI/GP is required. To set this, use: `default(parisize, 30000000000)`

Depending on the system specifications, it may be necessary to adjust the stage size or increase virtual memory (to simulate additional RAM). 

### Disk Space
Storage is required for saving and loading precomputed prime products. Each prime product of 10⁸ primes takes approximately 1GB. 

For example, `make_products(10, 8, 17)` requires 9.22GB of disk space. 

### Fermat primality test (Optional)
A program such as PrimeForm/GW can be useful for testing probable primes. 

Website: https://sourceforge.net/projects/openpfgw/

## Example Usage 
```
\\ Load the script.
read("sgsa1.0.gp");

\\ Set the PARI/GP available RAM to 28GB
default(parisize, 30000000000)

\\ Generate prime product stages of the first 10⁹ primes.
make_products(10, 8, 17)

\\ Save products to disk for faster future access.
save_products(17)

\\ Load saved products.
load_products(17)

\\ Generate prime candidates of 1,100,000 digits using an adapted Pocklington’s theorem and use SGSA batch mode of 17 stages.
pocklington_generator_batch("candidates", "2^3021377-1", 204068, 5, 6, 17, 1, 10000)

\\ \\ Performs the same operation using SGSA sequential mode.
pocklington_generator_seq("candidates", "2^3021377-1", 204068, 17, 1, 1000)

\\ Test the numbers in the file using SGSA batch mode.
file_batch("file.txt", 5, 6, 17)

\\ Test a single number with sequential mode.
sgsa(2^3021377-1,17)

\\ Factor a single number.
sgsa_factor(2^3021377-2,17)
```

## Commands
* `default(parisize, 30000000000);` Sets the PARI/GP available RAM in bytes (example: 30,000,000,000 bytes, approximately  28GB). 
* `time_update = 60000;` Sets the minimum interval (in milliseconds) between progress updates, which occur only when the candidate changes. (default: 60,000 ms).
* `log_info = 0;` Enables (1) or disables (0) detailed logging.
* `show_time = 0;` Enables (1) or disables (0) execution time display for each stage and candidate. 
* `calculate_digits = 1;` Enables (1) or disables (0) computing and displaying the number of digits in the prime and first candidate. (For large numbers, this may take time).
* `twinprime_search = 0;` Enables (1) or disables (0) twin prime search by testing both n and n–2 in the  pocklington_generator_batch(), pocklington_generator_seq(), and file_batch().

## Functions
### make_products(growth_factor, fixed_stage, max_stage)
**Description**: Computes and stores prime products across multiple stages.
* **Parameters**:
    * growth_factor: Integer > 1, determines the growth rate until the fixed_stage.
    * fixed_stage: Integer ≥ 1, specifies the number of the first fixed stages.
    * max_stage: Integer ≥ fixed_stage, sets the total of stages.
* **Output**: Displays prime product information and stores them in the internal vector stage_products.

The  growth_factor applies only up to the fixed_stage. From fixed_stage to max_stage, the stage size remains constant. 

As an example, `make_products(10, 3, 6)`. Results in:
* Stage 1: Primes 1 to 10.
* Stage 2: Primes 11 to 100.
* Stage 3: Primes 101 to 1,000.
* Stage 4: Primes 1,001 to 2,000.
* Stage 5: Primes 2,001 to 3,000.
* Stage 6: Primes 3,001 to 4,000.

Recommended usage:  `make_products(10,7,max_stage)` or `make_products(10,8,max_stage)`.

### save_products(max_stage)
**Description**: Saves the computed prime products vector stage_products to disk.
* **Parameters**:
    * max_stage: Integer ≥ 1, specifying how many stages to save.
* **Output**: Stores stage products in individual files named sgsa_stageprod_*.txt.
Ensure that any existing sgsa_stageprod_*.txt files are renamed or deleted before running this function. 

### load_products(max_stage)
**Description**: Loads prime product data from saved files into the internal vector stage_products. 
* **Parameters**:
    * max_stage: Integer ≥ 1, specifying how many stages to load.
* **Output**: Populates the internal vector stage_products with loaded data.

### pocklington_generator_batch("file_name", "q_string", e, batch_stage, divide_stage, max_stage, start_i, end_i)
**Description**: Generates adapted Pocklington's theorem prime candidates and uses SGSA batch mode.
* **Parameters**:
    * "file_name": String, base name for the output file, enclosed in quotes, e.g., "file".
    * "q_string": String, must be an expression representing a prime number, enclosed in quotes, e.g., "2^31-1".
    * e: Integer > 1, exponent for constructing large numbers.
    * batch_stage: Integer > 1, specifying the stage to begin batch mode.
    * divide_stage: 0 or Integer > batch_stage, specifying the stage at which division is used instead of multiplication. 0 disables this feature.
    * max_stage: Integer, specifying the maximum stage to apply.
    * start_i, end_i: Integers defining the iteration range for the batch.
* **Output**: Displays sieving information and saves prime candidates to a file.

The saved file can be used in PrimeForm/GW or other applications to perform a Fermat primality test and determine whether the candidate is a probable prime. 

To test with PrimeForm/GW: `pfgw64 file.txt`

This method applies an adapted Pocklington’s theorem to generate candidates using the following expression: 

`q * 10^e - 2 * i * q + 1` which is equivalent to `q * (10^e − 2 * i) + 1`.

If a probable prime is found, its primality can be verified efficiently using PrimeForm/GW.  To prove primality, create a helper file named helper.txt containing the q_string without quotes, which has to be prime. Then execute: 
```
pfgw64 -t -V -hhelper.txt -q"expression_for_n"
```
For twin primes, see the corresponding section. 

### pocklington_generator_seq("file_name", "q_string", e, max_stage, start_i, end_i)
**Description**: Generates Pocklington's theorem prime candidates and use SGSA sequential mode.
* **Parameters**:
    * "file_name": String, base name for the output file, enclosed in quotes, e.g., "file".
    * "q_string": String, must be an expression representing a prime number, enclosed in quotes, e.g., "prime".
    * e: Integer > 1, exponent for constructing large numbers.
    * max_stage: Integer, specifying the maximum stage to apply.
    * start_i, end_i: Integers defining the iteration range.
* **Output**: Displays sieving information and saves prime candidates to a file.

The saved file can be used in the same way as the files generated with pocklington_generator_batch().

For twin primes, see the corresponding section.

### file_batch("file_name",batch_stage,divide_stage,max_stage)
**Description**: Test the numbers in the file using SGSA batch mode.
* **Parameters**:
    * "file_name": String, name of the input file, enclosed in quotes, e.g., "file.txt".
    * batch_stage: Integer > 1, specifying the stage to begin batch mode.
    * divide_stage: 0 or Integer > batch_stage, specifying the stage at which division is used instead of multiplication. 0 disables this feature.
    * max_stage: Integer, specifying the maximum stage to apply.
* **Output**: Displays sieving information and saves prime candidates to a new file named by appending `_finished.txt` to the original file name.

### sgsa(n, max_stage)
**Description**: Test a single number with sequential mode.
* **Parameters**:
    * n: Integer or expression to be tested.
    * max_stage: Integer, specifying the maximum stage to apply.
* **Output**: Returns 1 if n passes all stages without finding a factor; otherwise, returns 0. 

### sgsa_factor(n, max_stage)
**Description**: Search for factors of a single number. 
* **Parameters**:
    * n: Integer or expression to be tested.
    * max_stage: Integer, specifying the maximum stage to apply.
* **Output**: Returns the found factors, if any.


## Twin primes

With pocklington_generator_batch(), pocklington_generator_seq() and file_batch(), it is possible to search for twin primes with this command `twinprime_search = 1`

This method applies SGSA to each candidate n and n-2. If either fails any sieve stage, the candidate is discarded. In this way, all candidates that pass all the stages will have done so for n and n-2. This significantly reduces the number of candidates that need to undergo the Fermat primality test.

If a probable twin prime is found, the primality of n-2 can be verified efficiently. To do this using PrimeForm/GW, create a helper file named `helper.txt` containing the q_string without quotes, which has to be prime. Then execute: 
```
pfgw64 -tp -V -hhelper.txt -q"n expression".
```
Make sure to use `-tp` instead of `-t` for n-2. 

## Extending the Code 
With knowledge of PARI/GP scripting, additional functions can be created to generate candidates using different methods, which can then be tested with SGSA via file_batch().

It may also be possible to adapt other sieving methods by incorporating parts of SGSA to improve performance. 

Pari/GP is particularly useful for handling extremely large variables, and its scripting language is easy to use. If the SGSA algorithm is properly integrated into a low-level language, it would likely be much more efficient.

## License

Copyright (c) 2025 Stage GCD Sieving Algorithm (SGSA)

This project is licensed under the Creative Commons Attribution 4.0 International License (CC BY 4.0). You are free to use, modify, and distribute the code, even commercially, as long as you give appropriate credit by mentioning the name of this script and linking to its source repository:

https://github.com/rrodriguezcunillera/sgsa

For the full license text, see: https://creativecommons.org/licenses/by/4.0/
