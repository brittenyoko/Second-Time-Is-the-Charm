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
    int m1[][] = {{3,2,4}, {5,1,9}, {2,3,0}};
    int m2[][] = {{3,8,7}, {5,6,2}, {1,7,9}};
    printMatrix(matrixMultiply(m1, m2)); // should be 23 64 61 / 29 109 118 / 21 34 20
  }
}
