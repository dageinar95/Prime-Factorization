# Prime Factorization
 
The program contains both a sequential and parallel algorithm for finding all prime numbers up to n (using the Sieve of Eratosthenes), and then factorizes the 100 largest primes in n^2. 

To run the program, first compile it using:

javac PrimeFactorization.java

Then run it like this:

java PrimeFactorization <n> <k>
 
where the Sieve of Eratosthenes will find all primes <= n and the program will factorize the 100 largest numbers in n*n. k is the number of threads to be used, and if k = 0 the program will use all available cores on the machine. 

An example of running the program:

java PrimeFactorization 2000000 0

The program will then find all primes <= 2 000 000 and factorize the 100 largest numbers of 4 000 000 000 000, using all available cores.
