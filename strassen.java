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
    int log2n = Math.log(n) / Math.log(2);
    if (floor(log2n) != ceil(log2n)) {
      Arrays.copyOf(m1, Math.pow(2, ceil(log2n)));
      Arrays.copyOf(m2, Math.pow(2, ceil(log2n)));
    }
    printMatrix(m1);
    printMatrix(m2);
    // int res[][] = new int[n][n];
    // for (int i = 0; i < n; i++) {
    //   for (int j = 0; j < n; j++) {
    //     res[i][j] = 0;
    //     for (int k = 0; k < n; k++) {
    //       res[i][j] += m1[i][k] * m2[k][j];
    //     }
    //   }
    // }
    // return res;
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
  }
}
