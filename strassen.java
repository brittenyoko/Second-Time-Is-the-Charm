package prog2;

import java.util.*;

public class Strassen {
  private static int[][] matrixMultiply(int[][] m1, int[][] m2) {
    int n = m1.length;
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

  private static int[][] strassenMultiply(int[][] m1, int[][] m2) {
    int n = m1.length;
    if (n == 1) {
      return {{m1[1][1] * m2[1][1]}};
    }
    double log2n = Math.log(n) / Math.log(2);
    if (Math.floor(log2n) != Math.ceil(log2n)) {
      n = (int) Math.pow(2, Math.ceil(log2n));
      int[][] m1_new = new int[n][n];
      int[][] m2_new = new int[n][n];
      for (int i = 0; i < m1.length; i++) {
        m1_new[i] = Arrays.copyOf(m1[i], n);
        m2_new[i] = Arrays.copyOf(m1[i], n);
      }
      m1 = m1_new;
      m2 = m2_new;
    }
    int[][][]
    return res;
  }

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

  private static void printMatrix(int[][] m) {
    int n = m.length;
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        System.out.print(m[i][j] + " ");
      }
      System.out.println();
    }
  }

  public static void main(String[] args) {
    int m[][] = {{3,2,4,6}, {5,1,9,4}, {2,3,0,1}, {2,1,8,4}};
    printMatrix(m);
    System.out.println();
    int[][][] submatrices = getSubmatrices(m);
    for (int i = 0; i < submatrices.length; i++) {
      printMatrix(submatrices[i]);
      System.out.println();
    }

    int m1[][] = {{3,2,4}, {5,1,9}, {2,3,0}};
    int m2[][] = {{1,0,2}, {6,7,1}, {3,9,0}};
    printMatrix(matrixMultiply(m1, m2));
    System.out.println();
    printMatrix(strassenMultiply(m1, m2));
  }
}
