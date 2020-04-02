package prog2;

import java.util.*;
import java.io.*;

public class Strassen {
  private static final SplittableRandom RAND = new SplittableRandom();
  private static final List<Integer> TEST_NS =
    new ArrayList<>(Arrays.asList(2,4,8,16,32,64,128,256,512,1024,2048));
  private static final int CROSSOVER = 136;
  private static final int TRIANGLES = 3;
  private static final int NORMAL_FLAG = 0;
  private static final int FIND_CROSSOVERS_FLAG = 1;
  private static final int COMPARE_STD_AND_STRASSEN_FLAG = 2;

  public static void main(String[] args) {
    // check args
    if (args.length == 0) {
      throw new IllegalArgumentException("Missing strassen argument(s)");
    }
    int flag = Integer.valueOf(args[0]);
    if (flag < 0 || flag > 5) {
      throw new IllegalArgumentException("Invalid flag");
    }
    if (flag == FIND_CROSSOVERS_FLAG) {
      findBestCrossovers();
    }
      else if (flag == COMPARE_STD_AND_STRASSEN_FLAG) {
      compareStrassenAndStd();
    }
    if (flag == TRIANGLES) {
        result();
    }
     else {
      if (args.length != 3) {
        throw new IllegalArgumentException("Missing or extra strassen arguments");
      }
      int dim = Integer.valueOf(args[1]);
      String fileName = args[2];
      if (dim < 0) {
        throw new IllegalArgumentException("Invalid dim arg");
      } else {
        // do strassen multiply on the matrices in the file and print the diagonal
        int[][][] matrices = txtFileToMatrices(fileName, dim);
        printDiagonal(strassenMultiply(matrices[0], matrices[1]));
      }
    }
  }

  // reads txt file into two matrices
  private static int[][][] txtFileToMatrices (String fileName, int dim) {
    int n_lines = 2 * (int) Math.pow(dim, 2);
    int half_lines = n_lines / 2;
    int[][][] matrices = new int[2][dim][dim];
    int[][] vals = new int[2][half_lines];
    try {
      // read file into matrices
      FileInputStream fstream = new FileInputStream(fileName);
      BufferedReader br = new BufferedReader(new InputStreamReader(fstream));
      int i = 0;
      String line = null;
      int matrix = 0;
      while (i < n_lines) {
        line = br.readLine();
        if (i < half_lines) {
          vals[0][i] = Integer.valueOf(line);
        } else {
          vals[1][i - half_lines] = Integer.valueOf(line);
        }
        i++;
      }
      fstream.close();
      matrices[0] = listToMatrix(vals[0]);
      matrices[1] = listToMatrix(vals[1]);
    } catch(Exception e) {
      e.printStackTrace();
    }
    return matrices;
  }

  // converts 1D list of matrix values to 2D matrix
  private static int[][] listToMatrix(int[] vals) {
    int n = (int) Math.sqrt(vals.length);
    int[][] m = new int[n][n];
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        m[i][j] = vals[i * n + j];
      }
    }
    return m;
  }

  // prints diagonal of matrix with trailing new line
  private static void printDiagonal(int[][] m) {
    for (int i = 0; i < m.length; i++) {
      System.out.println(m[i][i]);
    }
    System.out.println();
  }

  // prints a matrix
  private static void printMatrix(int[][] m) {
    int n = m.length;
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        System.out.print(m[i][j] + " ");
      }
      System.out.println();
    }
    System.out.println();
  }

  // standard matrix multiply
  private static int[][] stdMultiply(int[][] m1, int[][] m2) {
    int n = m1.length;
    if (n == 0) {
      return m1;
    }
    int res[][] = new int[n][n];
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        res[i][j] = 0;
        for (int k = 0; k < n; k++) {
          res[i][j] += m1[i][k] * m2[k][j];
        }
      }
    }
    return res;
  }

  // strassen multiplication with default crossover point
  private static int[][] strassenMultiply(int[][] m1, int[][] m2) {
    return strassenMultiply(m1, m2, CROSSOVER);
  }

  // strassen multiplication with given crossover point
  private static int[][] strassenMultiply(int[][] m1, int[][] m2, int crossoverPoint) {
    int n = m1.length;
    int orig_n = n;
    // handle base cases
    if (n == 0) {
      return m1;
    }
    if (n == 1) {
      int res[][] = {{m1[0][0] * m2[0][0]}};
      return res;
    }
    if (n <= crossoverPoint) {
      return stdMultiply(m1, m2);
    }
    // pad matrices with zeros if necessary to reach even size
    boolean needToPad = n % 2 != 0;
    if (needToPad) {
      n += 1;
      m1 = padMatrix(m1, n);
      m2 = padMatrix(m2, n);
    }
    int[][][] subs1 = getSubmatrices(m1);
    int[][][] subs2 = getSubmatrices(m2);
    int[][][] prods = new int[7][n][n];
    // get products
    int pt = crossoverPoint;
    prods[0] = strassenMultiply(subs1[0], sub(subs2[1], subs2[3]), pt);
    prods[1] = strassenMultiply(add(subs1[0], subs1[1]), subs2[3], pt);
    prods[2] = strassenMultiply(add(subs1[2], subs1[3]), subs2[0], pt);
    prods[3] = strassenMultiply(subs1[3], sub(subs2[2], subs2[0]), pt);
    prods[4] = strassenMultiply(add(subs1[0], subs1[3]), add(subs2[0], subs2[3]), pt);
    prods[5] = strassenMultiply(sub(subs1[1], subs1[3]), add(subs2[2], subs2[3]), pt);
    prods[6] = strassenMultiply(sub(subs1[0], subs1[2]), add(subs2[0], subs2[1]), pt);
    // get resulting matrix
    int[][][] res_subs = new int[4][n / 2][n / 2];
    res_subs[0] = add(sub(add(prods[4], prods[3]), prods[1]), prods[5]);
    res_subs[1] = add(prods[0], prods[1]);
    res_subs[2] = add(prods[2], prods[3]);
    res_subs[3] = sub(sub(add(prods[0], prods[4]), prods[2]), prods[6]);
    int[][] res = combineSubmatrices(res_subs);
    return needToPad ? truncateMatrix(res, orig_n) : res;
  }

  // adds 2 matrices
  private static int[][] add(int[][] m1, int[][] m2) {
    return addOrSub(false, m1, m2);
  }

  // subtracts 2 matrices
  private static int[][] sub(int[][] m1, int[][] m2) {
    return addOrSub(true, m1, m2);
  }

  // adds or subtracts 2 matrices based on given flag
  private static int[][] addOrSub(boolean subtract, int[][] m1, int[][] m2) {
    int n = m1.length;
    int res[][] = new int[n][n];
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        res[i][j] = m1[i][j] + (subtract ? -1 : 1) * m2[i][j];
      }
    }
    return res;
  }

  // pads matrix with zeros until it reaches size n
  private static int[][] padMatrix(int[][] m, int n) {
    int[][] m_new = new int[n][n];
    for (int i = 0; i < m.length; i++) {
      m_new[i] = Arrays.copyOf(m[i], n);
    }
    return m_new;
  }

  // truncates matrix until it reaches specified size
  private static int[][] truncateMatrix(int[][] m, int n) {
    int[][] m_new = new int[n][n];
    for (int i = 0; i < n; i++) {
      m_new[i] = Arrays.copyOf(m[i], n);
    }
    return m_new;
  }

  // gets the 4 submatrices for Strassen's algorithm
  private static int[][][] getSubmatrices(int[][] m) {
    int n = m.length;
    int n1 = n / 2;
    int[][][] submatrices = new int[4][n1][n1];
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        if (i < n1) {
          if (j < n1) {
            submatrices[0][i][j] = m[i][j];
          } else {
            submatrices[1][i][j - n1] = m[i][j];
          }
        } else {
          if (j < n1) {
            submatrices[2][i - n1][j] = m[i][j];
          } else {
            submatrices[3][i - n1][j - n1] = m[i][j];
          }
        }
      }
    }
    return submatrices;
  }

  // combines 4 submatrices into 1
  private static int[][] combineSubmatrices(int[][][] submatrices) {
    int n1 = submatrices[0].length;
    int n = n1 * 2;
    int[][] m = new int[n][n];
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        if (i < n1) {
          if (j < n1) {
            m[i][j] = submatrices[0][i][j];
          } else {
            m[i][j] = submatrices[1][i][j - n1];
          }
        } else {
          if (j < n1) {
            m[i][j] = submatrices[2][i - n1][j];
          } else {
            m[i][j] = submatrices[3][i - n1][j - n1];
          }
        }
      }
    }
    return m;
  }

  // gets matrix of size n randomly filled with 0s or 1s
  private static int[][] getRandomMatrix(int n) {
    int[][] m = new int[n][n];
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        m[i][j] = RAND.nextInt(10) % 2;
      }
    }
    return m;
  }

  // run strassens on different sized matrices with different crossover pts to find the optimal
  private static void findBestCrossovers() {
    Map<Integer, List<Integer>> best_crossovers = new LinkedHashMap<>();
    int trials = 15;
    int incr = 4;
    int percent_within = 5;
    for (int n : TEST_NS) {
      // create random matrices to test on
      int[][][] matrices = new int[trials * 2][n][n];
      for (int i = 0; i < matrices.length; i++) {
        matrices[i] = getRandomMatrix(n);
      }
      List<Double> avg_times = new ArrayList<>();
      double min_avg_time = Double.MAX_VALUE;
      int start = n >= 128 ? (n >= 1024 ? 100 : 40) : 0;
      int end = n >= 1024 ? 300 : 256;
      // test different crossover points in increments of 'incr'
      for (int crossover = start; crossover < Math.min(end, n + incr); crossover += incr) {
        System.out.println(String.format("TESTING %dx%d MATRIX, CROSSOVER AT SIZE %d...",
          n, n, crossover));
        List<Double> times = new ArrayList<>();
        for (int i = 0; i < trials; i++) {
          System.out.print(".");
          double startTime = System.nanoTime();
          strassenMultiply(matrices[i * 2], matrices[i * 2 + 1], crossover);
          times.add((System.nanoTime() - startTime) / 1e6);
        }
        // get average time of all the trials
        double avg = times.stream().mapToDouble(a -> a).average().getAsDouble();
        System.out.println(String.format("AVG TIME: %.03f ms\n", avg));
        avg_times.add(avg);
        min_avg_time = Math.min(min_avg_time, avg);
      }
      // find the crossover point that resulted in the lowest avg time
      int min_index = avg_times.indexOf(min_avg_time);
      int best_crossover = start + incr * min_index;
      best_crossover = best_crossover > n ? n : best_crossover;
      // calculate range of crossovers that result in time within certain percent of best
      List<Integer> crossovers_within = new ArrayList<>();
      for (int i = 0; i < avg_times.size(); i++) {
        double time = avg_times.get(i);
        if ((time - min_avg_time) / min_avg_time < percent_within * .01) {
          crossovers_within.add(start + incr * i);
        }
      }
      Collections.sort(crossovers_within);
      crossovers_within.add(0, best_crossover);
      best_crossovers.put(n, crossovers_within);
      System.out.println();
    }
    System.out.println("BEST CROSSOVERS:\n" + crossoversMapToStr(best_crossovers));
  }

  // use both strassens and standard mult. on different size matrices and see which is faster
  private static void compareStrassenAndStd() {
    int trials = 10;
    for (int n : TEST_NS) {
      System.out.print(String.format("SIZE %d MATRIX", n));
      // create random matrices to test on
      int[][][] matrices = new int[trials * 2][n][n];
      for (int i = 0; i < matrices.length; i++) {
        matrices[i] = getRandomMatrix(n);
      }
      List<Double> std_times_ms = new ArrayList<>();
      List<Double> strassen_times_ms = new ArrayList<>();
      for (int i = 0; i < trials; i++) {
        // perform standard and strassen multiplication on the matrices
        System.out.print(".");
        double stdMultTimeMs;
        if (n >= 2048) {
          if (i == 0) {
            std_times_ms.add(Double.MAX_VALUE);
          }
        } else {
          double startTime1 = System.nanoTime();
          stdMultiply(matrices[i * 2], matrices[i * 2 + 1]);
          std_times_ms.add((System.nanoTime() - startTime1) / 1e6);
        }
        double startTime2 = System.nanoTime();
        strassenMultiply(matrices[i * 2], matrices[i * 2 + 1]);
        strassen_times_ms.add((System.nanoTime() - startTime2) / 1e6);
      }
      // get average time of the trials for both methods and compare the two
      double std_avg = std_times_ms.stream().mapToDouble(a -> a).average().getAsDouble();
      double strassen_avg = strassen_times_ms.stream().mapToDouble(a -> a).average().getAsDouble();
      System.out.println(String.format("\nTRADITIONAL AVG TIME: %s",
        std_avg == Double.MAX_VALUE ? "N/A" : String.format("%.03f ms", std_avg)));
      System.out.println(String.format("STRASSEN AVG TIME: %.03f ms", strassen_avg));
      if (Math.abs((std_avg - strassen_avg) / std_avg) < 0.1) {
        System.out.println("  ABOUT EQUAL\n");
      } else {
        System.out.println(String.format("  %s faster\n",
          std_avg < strassen_avg ? "TRADITIONAL" : "STRASSEN"));
      }
    }
  }

  // converts map from matrix size to best range of crossovers to a string for printing
  private static String crossoversMapToStr(Map<Integer, List<Integer>> map) {
    String str = "";
    for (Map.Entry<Integer, List<Integer>> entry : map.entrySet()) {
      int n = entry.getKey();
      List<Integer> crossovers = entry.getValue();
      int best_crossover = crossovers.get(0);
      crossovers.remove(0);
      str += String.format("%d -> %d, %s\n", n, best_crossover, crossovers.toString());
    }
    return str;
  }

  // checks that two matrices are equal
  private static boolean areEqual(int[][] m1, int[][] m2) {
    for (int i = 0; i < m1.length; i++) {
      if (!Arrays.equals(m1[i], m2[i])) {
        return false;
      }
    }
    return true;
  }

  // Creates random graph of 1024 vertices with some given probability p
  private static int[][] graphMaker(double p){
    int[][] mat_rand = new int[1024][1024];
    for (int i = 0; i < 1024; i++) {
      for (int j = 0; j < 1024; j++) {
        int rand = RAND.nextInt(100);
        mat_rand[i][j] = rand < p * 100 ? 1 : 0;
      }
    }
    return mat_rand;
  }

  // Takes in a random graph of 1024 vertices and returns the number of triangles
  private static float triangles(int[][] mat_rand){
      int[][] A = mat_rand;
      int[][] B = strassenMultiply(A,A);
      int[][] C = strassenMultiply(A,B);
      float sum = 0;
      for (int i = 0; i < 1024; i++) {
        for (int j = 0; j < 1024; j++) {
          if (i == j) {
            sum = sum + C[i][j];
          }
        }
      }
      return sum / 6;
  }
  private static void result() {
    for (int i = 0; i < 5; i++) {
      Float flt = triangles(graphMaker(.03));
      System.out.println(flt);
    }
  }
}
