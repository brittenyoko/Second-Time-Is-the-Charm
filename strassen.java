package prog2;

import java.util.*;

public class Strassen {
  public static int crossoverPoint;
  public static final int[] TEST_NS = {8,16,32,64,128,256};

  // traditional matrix multiply
  private static int[][] matrixMultiply(int[][] m1, int[][] m2) {
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

  // strassen multiplication
  private static int[][] strassenMultiply(int[][] m1, int[][] m2) {
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
      return matrixMultiply(m1, m2);
    }
    // pad matrices with zeros if necessary to reach size that is power of 2
    double log2n = Math.log(n) / Math.log(2);
    boolean needToPad = Math.floor(log2n) != Math.ceil(log2n);
    if (needToPad) {
      n = (int) Math.pow(2, Math.ceil(log2n));
      m1 = padMatrix(m1, n);
      m2 = padMatrix(m2, n);
    }
    int[][][] subs1 = getSubmatrices(m1);
    int[][][] subs2 = getSubmatrices(m2);
    int[][][] prods = new int[7][n][n];
    // get products
    prods[0] = strassenMultiply(subs1[0], sub(subs2[1], subs2[3]));
    prods[1] = strassenMultiply(add(subs1[0], subs1[1]), subs2[3]);
    prods[2] = strassenMultiply(add(subs1[2], subs1[3]), subs2[0]);
    prods[3] = strassenMultiply(subs1[3], sub(subs2[2], subs2[0]));
    prods[4] = strassenMultiply(add(subs1[0], subs1[3]), add(subs2[0], subs2[3]));
    prods[5] = strassenMultiply(sub(subs1[1], subs1[3]), add(subs2[2], subs2[3]));
    prods[6] = strassenMultiply(sub(subs1[0], subs1[2]), add(subs2[0], subs2[1]));
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

  // pads matrix with zeros until it reaches specified size
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

  // gets matrix of size n randomly filled with 0s or 1s
  public static int[][] getRandomMatrix(int n) {
    SplittableRandom rand = new SplittableRandom();
    int[][] m = new int[n][n];
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        m[i][j] = rand.nextInt(10) % 2;
      }
    }
    return m;
  }

  public static void findCrossoverPoint() {
    for (int i = 0; i < TEST_NS.length; i++) {
      int n = TEST_NS[i];
      System.out.println(String.format("SIZE %d MATRIX", n));
      int[][] m1 = getRandomMatrix(n);
      int[][] m2 = getRandomMatrix(n);
      // time traditional and strassen multiplication and compare the two
      double startTime1 = System.nanoTime();
      int[][] normal_res = matrixMultiply(m1, m2);
      double normalMultTimeMs = (System.nanoTime() - startTime1) / 1e6;

      double startTime2 = System.nanoTime();
      int[][] strassen_res = strassenMultiply(m1, m2);
      double strassenTimeMs = (System.nanoTime() - startTime2) / 1e6;

      System.out.println(String.format("Normal time: %.06f ms", normalMultTimeMs));
      System.out.println(String.format("Strassen time: %.06f ms", strassenTimeMs));
      System.out.println(String.format("%s faster",
        normalMultTimeMs < strassenTimeMs ? "NORMAL" : "STRASSEN"));
      System.out.println();
    }
  }

  public static void main(String[] args) {
    // int m1[][] = {{3,2,4},{5,1,9},{2,3,0}};
    // int m2[][] = {{1,0,2},{6,7,1},{3,9,0}};
    // printMatrix(matrixMultiply(m1, m2));
    // printMatrix(strassenMultiply(m1, m2));
    findCrossoverPoint();
  }
}
