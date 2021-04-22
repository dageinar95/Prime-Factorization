import java.util.*;
import java.util.concurrent.locks.ReentrantLock;
import java.util.concurrent.locks.Lock;
import java.math.BigDecimal;
import java.math.RoundingMode;

class PrimeFactorization {

  class SequentialSOE { //this class includes the sequential factorization

    /**
     * Declaring all the global variables
     *
     */
    int n, root, numOfPrimes;
    byte[] oddNumbers;

    /**
     * Constructor that initializes the global variables
     * @param n Prime numbers up until (and including if prime) 'n' is found
     */
    SequentialSOE(int n) {
      this.n = n;
      root = (int) Math.sqrt(n);
      oddNumbers = new byte[(n / 16) + 1];
    }

    /**
     * Performs the sieve and collects the primes produced by the sieve.
     * @return An array containing all the primes up to and including 'n'.
     */
    int[] getPrimes() {
      if (n <= 1) return new int[0];

      sieve();

      return collectPrimes();
    }

    /**
     * Iterates through the array to count the number of primes found,
     * creates an array of that size and populates the new array with the primes.
     * @return An array containing all the primes up to and including 'n'.
     */
    private int[] collectPrimes() {

      int start = (root % 2 == 0) ? root + 1 : root + 2; //starting at the root of n, since we have already counted the primes up to the root when sieving

      for (int i = start; i <= n; i += 2)
        if (isPrime(i))
          numOfPrimes++;

      int[] primes = new int[numOfPrimes];

      primes[0] = 2;

      int j = 1;

      for (int i = 3; i <= n; i += 2)
        if (isPrime(i))
          primes[j++] = i;

      return primes;
    }

    /**
     * Performs the Sieve Of Eratosthenes
     */
    private void sieve() {
      mark(1);
      numOfPrimes = 1;
      int prime = nextPrime(1);

      while (prime != -1) {
        traverse(prime);
        prime = nextPrime(prime);
        numOfPrimes++;
      }
    }

    /**
     * Marks all odd number multiples of 'prime', starting from prime * prime.
     * @param prime The prime used to mark the composite numbers.
     */
    private void traverse(int prime) {
      for (int i = prime*prime; i <= n; i += prime * 2)
        mark(i);
    }

    /**
     * Finds the next prime in the sequence. If there are no more left, it
     * simply returns -1.
     * @param  prev The last prime that has been used to mark all non-primes.
     * @return      The next prime or -1 if there are no more primes.
     */
    private int nextPrime(int prev) {
      for (int i = prev + 2; i <= root; i += 2)
        if (isPrime(i))
          return i;

      return -1;
    }

    /**
     * Checks if a number is a prime number. If 'num' is prime, it returns true.
     * If 'num' is composite, it returns false.
     * @param  num The number to check.
     * @return     A boolean; true if prime, false if not.
     */
    private boolean isPrime(int num) {
      int bitIndex = (num % 16) / 2;
      int byteIndex = num / 16;

      return (oddNumbers[byteIndex] & (1 << bitIndex)) == 0;
    }

    /**
     * Marks the number 'num' as a composite number (non-prime)
     * @param num The number to be marked non-prime.
     */
    private void mark(int num) {
      int bitIndex = (num % 16) / 2;
      int byteIndex = num / 16;
      oddNumbers[byteIndex] |= (1 << bitIndex);
    }

    /**
     * Prints the primes found.
     * @param primes The array containing all the primes.
     */
    private void printPrimes(int[] primes) {
      for (int prime : primes)
        System.out.println(prime);
    }

    private void printFirstXPrimes(int[] primes, int X) { //method for printing the primes that are less than X
      System.out.println("All primes less than " + X + ": ");
      int counter = 0;
      while (true) {
        if (counter < primes.length && primes[counter] < X) {
          System.out.print(primes[counter] + ", ");
          counter++;
        }
        else {
          break;
        }
      }
      System.out.println();
    }

    // method for factorizing all largest numbers
    // private void factorize(int[] primes) {
    //   int[] largest = new int[100]; //array to contain the largest numbers in n*n
    //
    //   int index = 0;
    //   for (int i = (n*n) - 99; i <= n*n; i++) { //saving the 100 largest numbers to the array
    //     largest[index] = i;
    //     index++;
    //   }
    //
    //   for (int i = 99; i < largest.length; i++) {
    //     int temp = largest[i];
    //     int counter = 0;
    //     while (counter < primes.length) {
    //       if (temp % primes[counter] == 0) { //if divisible by a prime, starting at the smallest one
    //         temp = temp / primes[counter];
    //         System.out.println("Factor for " + largest[i] + ": " + primes[counter]);
    //       }
    //       else { //if not divisible, increase counter and try next prime
    //         counter++;
    //       }
    //     }
    //   }
    // }

    // function for factorizing a single number and returning an arraylist of factors used
    private ArrayList<Integer> factorize(long num, int[] primes) {
      ArrayList<Integer> factors = new ArrayList<Integer>();
      long temp = num;

      int counter = 0;
      while (counter < primes.length) {
        if (temp % primes[counter] == 0) { //if divisible by a prime, starting at the smallest one
          temp = temp / primes[counter];
          factors.add(primes[counter]);
        }
        else { //if not divisible, increase counter and try next prime
          counter++;
        }
      }

      return factors;
    }

  }


  class ParallelSOE {

    int n, k, root, numOfPrimes;
    byte[] oddNumbers;

    ParallelSOE(int n, int k) {
      this.n = n;
      this.k = k;
      root = (int) Math.sqrt(n);
      oddNumbers = new byte[(n/16) + 1];
    }

    class Worker1 implements Runnable { //Worker1 will handle the parallelization of the sieve
      int id, start, end;
      ReentrantLock rl = new ReentrantLock();

      Worker1(int id, int start, int end) {
        this.id = id;
        this.start = start; //which byte index to start at
        this.end = end; //byte index to end before
      }

      public void run() {
        int root = (int) Math.sqrt(end*16);
        SequentialSOE soe = new PrimeFactorization().new SequentialSOE(root);
        int[] primes = soe.getPrimes(); //generating all the primes less than the root of the highest number the thread has to cover

        int low = start*16;
        int high = (end*16)-1;
        int temp = low;
        // System.out.println("Thread " + id + ": low is " + low + ", high is " + high);

        if (id == 0) { //for the first thread we already know the "pattern"
          mark(1);
          temp = nextPrime(1, high);

          while (temp <= root && temp != -1) {
            traverse(temp, low, high);
            temp = nextPrime(temp, high);
          }
        }

        else {
          for (int prime : primes) {
            if (prime != 2) {
              traverse(prime, low, high);
            }
          }
        }

        //updating the numOfPrimes variable
        if (id == k-1) { //the last thread will use all the primes up to the root of n
          numOfPrimes += primes.length;
        }
      }


    }

    private void traverse(int prime, int start, int end) { //modified traverse, to begin marking of non-primes at the correct index.
      int newStart = prime * prime;
      while (newStart < start) {
        newStart += prime * 2;
      }
      for (int i = newStart; i <= end; i += prime * 2) {
        mark(i);
      }
    }

    private int[] getPrimes() {
      int nums = n / k; //how many numbers each thread will cover
      int bytes = oddNumbers.length / k; //how many bytes each thread, except the last, will cover
      Thread[] threads = new Thread[k];

      for (int i = 0; i < k-1; i++) {
        threads[i] = new Thread(new Worker1(i, i*bytes, (i+1)*bytes));
      }

      threads[k-1] = new Thread(new Worker1(k-1, (k-1)*bytes, oddNumbers.length)); //last thread takes the "rest" - useful if the amount of bytes doesn't add up

      numOfPrimes = 1; //since each thread will skip even numbers, start the variable at 1 to account for 2

      for (Thread t : threads) { t.start(); }
      try {
        for (Thread t : threads) { t.join(); }
      } catch (Exception e) {
        e.printStackTrace();
      }


      return collectPrimes();
    }

    private int[] collectPrimes() {

      int start = (root % 2 == 0) ? root + 1 : root + 2; //starting at the root of n, since we have already counted the primes up to the root when sieving

      for (int i = start; i <= n; i += 2) //counting the rest of the primes that each thread found
        if (isPrime(i)) {
          numOfPrimes++;
        }

      int[] primes = new int[numOfPrimes-1];
      // System.out.println("numOfPrimes = " + primes.length);

      primes[0] = 2;

      int j = 1;

      for (int i = 3; i <= n; i += 2)
        if (isPrime(i)) {
          primes[j++] = i;
        }

      return primes;
    }

    private int nextPrime(int prev, int high) { //modified nextPrime, with highest number
      for (int i = prev + 2; i <= high; i += 2)
        if (isPrime(i))
          return i;

      return -1;
    }

    private boolean isPrime(int num) {
      int bitIndex = (num % 16) / 2;
      int byteIndex = num / 16;

      return (oddNumbers[byteIndex] & (1 << bitIndex)) == 0;
    }

    private void mark(int num) {
      int bitIndex = (num % 16) / 2;
      int byteIndex = num / 16;
      oddNumbers[byteIndex] |= (1 << bitIndex);
    }

    private void printFirstXPrimes(int[] primes, int X) { //method for printing the primes that are less than X
      System.out.println("All primes less than " + X + ": ");
      int counter = 0;
      while (true) {
        if (counter < primes.length && primes[counter] < X) {
          System.out.print(primes[counter] + ", ");
          counter++;
        }
        else {
          break;
        }
      }
      System.out.println();
    }

  }

  class ParallelFac {
    int[] primes;
    int n, k;

    ParallelFac(int[] primes, int n, int k) {
      this.primes = primes; //variable to store all primes generated from the sieve
      this.n = n;
      this.k = k;
    }

    class Worker2 implements Runnable { //Worker2 will handle the parallelization of the factorization
      int id, firstPrime;
      long num;
      Lock rl;
      ArrayList<Integer> factors;

      Worker2(int id, int firstPrime, long num, Lock rl, ArrayList<Integer> factors) {
        this.id = id;
        this.firstPrime = firstPrime; //index of the first prime the thread will use. the next index will be firstPrime += k
        this.num = num;
        this.rl = rl;
        this.factors = factors;
      }

      public void run() {
        int nextPrime = firstPrime;

        while (nextPrime < primes.length) {
          if (num % primes[nextPrime] == 0) {
            try {
              rl.lock();
              num = num / primes[nextPrime];
              factors.add(primes[nextPrime]);
            }
            catch (Exception e) {
              e.printStackTrace();
            }
            finally {
              rl.unlock();
            }
          }
          else {
            nextPrime += k;
          }
        }
      }
    }

    private ArrayList<Integer> factorize(long n) {
      ArrayList<Integer> factors = new ArrayList<Integer>();
      long num = n;
      int root = (int) Math.sqrt(n);
      int primesPerThread = primes.length / k;
      Lock rl = new ReentrantLock();

      Thread[] threads = new Thread[k];

      for (int i = 0; i < k; i++) {
        threads[i] = new Thread(new Worker2(i, i, num, rl, factors));
      }

      for (Thread t : threads) { t.start(); }

      try {
        for (Thread t : threads) { t.join(); }
      } catch (Exception e) {
        e.printStackTrace();
      }

      return factors;
    }


  }

  public static void main(String[] args) {
    int n, k;

    try {
      n = Integer.parseInt(args[0]);
      k = Integer.parseInt(args[1]);
      if (n <= 16) throw new Exception();
      else if (k < 0) throw new Exception();
    } catch(Exception e) {
      System.out.println("Correct use of program is: " +
      "java PrimeFactorization <n> <k> where <n> is an integer > 16 and k is the number of threads to be used (0 if equal to cores available).");
      return;
    }

    if (k == 0) {
      System.out.println("k = 0 --> Using all available cores.");
      k = Runtime.getRuntime().availableProcessors();
    }

    long[] largest = new long[100]; //array to contain the 100 largest numbers in n*n
    int index = 0;
    long nSquared = (long) n * (long) n;
    for (long i = nSquared - 99; i <= nSquared; i++) { //saving the 100 largest numbers to the array
      largest[index] = i;
      index++;
    }

    Oblig3Precode precode = new Oblig3Precode(n);

    ////////////Sequential execution:

    System.out.println("\nBeginning sequential Sieve of Eratosthenes... ");
    SequentialSOE soe = new PrimeFactorization().new SequentialSOE(n); //Getting all the primes equal to and below 'n'

    int[] primes = null;
    double[] runtimes1 = new double[7];
    long start, end;
    start = end = 0;
    for (int i = 0; i < 7; i++) {
      start = System.nanoTime();
      primes = soe.getPrimes();
      end = System.nanoTime();
      runtimes1[i] = (end-start)/1000000.0;
      System.out.println("Run " + i + " sieve time: " + runtimes1[i]);
    }
    Arrays.sort(runtimes1);

    soe.printFirstXPrimes(primes, 100); //printing the primes < 100
    System.out.println("Sequential sieve of Eratosthenes median execution time: " + runtimes1[3] + " milliseconds.");

    System.out.println("\nBeginning sequential factorization... ");

    double[] runtimes2 = new double[7];
    long factorTime;
    for (int j = 0; j < 7; j++) { //7 runs
      factorTime = 0;
      for (int i = 0; i < largest.length; i++) { //iterating each of the 100 largest elements
        start = System.nanoTime();
        ArrayList<Integer> factors = soe.factorize(largest[i],primes);
        end = System.nanoTime();
        factorTime += (end-start); //adding only the time it took to factorize each number

        if (j == 0) { //we only need to save the factors to file once (one run)
          //not saving here, saving in the parallel exectution instead
          // saveFactors(factors, precode, largest[i], n);
        }
      }
      System.out.println("Run " + j + " factorization time: " + factorTime/1000000.0);
      runtimes2[j] = factorTime/1000000.0;
    }
    Arrays.sort(runtimes2);

    System.out.println("Sequential factorization median execution time: " + runtimes2[3] + " milliseconds.");

    ////////////Parallel execution:

    System.out.println("\nBeginning parallel Sieve of Eratosthenes... ");

    ParallelSOE psoe = new PrimeFactorization().new ParallelSOE(n,k);
    int[] primes2 = null;
    double[] runtimes3 = new double[7];
    for (int i = 0; i < 7; i++) {
      start = System.nanoTime();
      primes2 = psoe.getPrimes();
      end = System.nanoTime();
      runtimes3[i] = (end-start)/1000000.0;
      System.out.println("Run " + i + " parallel sieve time: " + runtimes3[i]);
      System.out.println("Run " + i + " speedup: " + round(runtimes1[3]/runtimes3[i],2)); //comparing the run's execution time to the sequential median
    }
    Arrays.sort(runtimes3);

    psoe.printFirstXPrimes(primes2, 100); //printing the primes < 100
    System.out.println("Parallel sieve of Eratosthenes median execution time: " + runtimes3[3] + " milliseconds.");

    if (Arrays.equals(primes, primes2)) {
      System.out.println("Success: Sequential and parallel sieves produced the same primes!");
    }
    else {
      System.out.println("Failure: Sequential and parallel sieves did not produce the same primes.");
    }

    System.out.println("\nBeginning parallel factorization... ");

    ParallelFac pf = new PrimeFactorization().new ParallelFac(primes2, n, k);

    double[] runtimes4 = new double[7];
    long factorTime2;
    for (int j = 0; j < 7; j++) { //7 runs
      factorTime2 = 0;
      for (int i = 0; i < largest.length; i++) { //iterating each of the 100 largest elements
        start = System.nanoTime();
        ArrayList<Integer> factors2 = pf.factorize(largest[i]);
        end = System.nanoTime();
        factorTime2 += (end-start); //adding only the time it took to factorize each number

        if (j == 0) { //we only need to save the factors to file once (one run)
          saveFactors(factors2, precode, largest[i], n);
        }
      }
      runtimes4[j] = factorTime2/1000000.0;
      System.out.println("Run " + j + " factorization time: " + runtimes4[j]);
      System.out.println("Run " + j + " speedup: " + round(runtimes2[3]/runtimes4[j],2));
    }
    Arrays.sort(runtimes4);

    System.out.println("Parallel factorization median execution time: " + runtimes4[3] + " milliseconds.");


  }

  public static void saveFactors(ArrayList<Integer> factors, Oblig3Precode precode, long num, int n) { //method for saving factors to file using the precode
    if (factors.isEmpty()) { //if there are no factors for a number, then it is a prime and is added as its own factor
      precode.addFactor(num, num);
    }
    else {
      long temp = 1; //variable to check that the factorization is complete
      for (int factor : factors) {
        temp = temp*factor; //"rebuilding" the number (largest[i])
        precode.addFactor(num, factor);
      }
      if (temp != num) { //if temp is different from the number, then the remaining prime is added as a factor
        precode.addFactor(num, (num/temp));
      }
    }
    precode.writeFactors();
  }

  public static double round(double value, int places) { //method for rounding double, used for comparison with SeqNotTransposed
    if (places < 0) throw new IllegalArgumentException();

    BigDecimal bd = BigDecimal.valueOf(value);
    bd = bd.setScale(places, RoundingMode.HALF_UP);
    return bd.doubleValue();
  }

}
