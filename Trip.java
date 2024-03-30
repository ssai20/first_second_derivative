import java.util.*;

public class Trip {


    public static int[][] mapAreaToMatrix(String area, String it)
    {
        int[][] matrix = new int[4][4];
        int[] simbolsToArray = new int[matrix.length*matrix.length];
        if (it.equals("Human"))
        {
            for (int i=0;i<simbolsToArray.length;i++)
            {
                if (area.charAt(i) == 'S')
                    simbolsToArray[i] = 5;
                if (area.charAt(i) == 'W')
                    simbolsToArray[i] = 2;
                if (area.charAt(i) == 'T')
                    simbolsToArray[i] = 3;
                if (area.charAt(i) == 'P')
                    simbolsToArray[i] = 1;
            }
        }
        if (it.equals("Swamper"))
        {
            for (int i=0;i<simbolsToArray.length;i++)
            {
                if (area.charAt(i) == 'S')
                    simbolsToArray[i] = 2;
                if (area.charAt(i) == 'W')
                    simbolsToArray[i] = 2;
                if (area.charAt(i) == 'T')
                    simbolsToArray[i] = 5;
                if (area.charAt(i) == 'P')
                    simbolsToArray[i] = 2;
            }
        }
        if (it.equals("Woodman"))
        {
            for (int i=0;i<simbolsToArray.length;i++)
            {
                if (area.charAt(i) == 'S')
                    simbolsToArray[i] = 3;
                if (area.charAt(i) == 'W')
                    simbolsToArray[i] = 3;
                if (area.charAt(i) == 'T')
                    simbolsToArray[i] = 2;
                if (area.charAt(i) == 'P')
                    simbolsToArray[i] = 2;
            }
        }

        int k=0;
        for (int i=0;i<matrix.length;i++)
        {
            for (int j=0;j<matrix.length;j++)
            {
                matrix[i][j] = simbolsToArray[k];
                k++;
            }
        }
        return matrix;
    }


    public static int searchAlgorithm(int[][] pathMatrix, int iStart, int jStart, int iEnd, int jEnd)
    {
        int[][] markers = new int[pathMatrix.length][pathMatrix.length];
        boolean[][] usedNode = new boolean[pathMatrix.length][pathMatrix.length];
        for (int i=0;i<markers.length;i++)
        {
            for (int j=0; j<markers.length;j++)
            {
                markers[i][j] = 1000;
                usedNode[i][j] = false;
            }
        }
        markers[iStart][jStart] = 0;


        for (int i=iStart;i<pathMatrix.length;i++)
        {
            for (int j=jStart;j<pathMatrix.length;j++)
            {
                if (j-1>=0) {
                    if ((!usedNode[i][j - 1]) && (markers[i][j - 1] > markers[i][j] + pathMatrix[i][j - 1])) {
                        markers[i][j - 1] = markers[i][j] + pathMatrix[i][j - 1];
                    }
                }
                if (i-1>=0)
                {
                    if ((!usedNode[i - 1][j]) && (markers[i - 1][j] > markers[i][j] + pathMatrix[i - 1][j])) {
                        markers[i - 1][j] = markers[i][j] + pathMatrix[i - 1][j];
                    }
                }

                if (j+1<pathMatrix.length)
                {
                    if ((!usedNode[i][j + 1]) && (markers[i][j+1] > markers[i][j]+pathMatrix[i][j+1]))
                        {
                            markers[i][j + 1] = markers[i][j] + pathMatrix[i][j + 1];
                        }
                }
                if (i+1<pathMatrix.length)
                {
                    if ((!usedNode[i + 1][j]) && (markers[i + 1][j] > markers[i][j] + pathMatrix[i + 1][j])) {
                        markers[i + 1][j] = markers[i][j] + pathMatrix[i + 1][j];
                    }
                }
                usedNode[i][j] = true;
            }
        }

        return markers[iEnd][jEnd];
    }


    public static int findingMinimalPath(String area, String it)
    {
        if ( (area.length()>16) || (area.length()<=0) ) throw new RuntimeException("Path should have length equals 16");
        int[][] pathMatrix = mapAreaToMatrix(area, it);
        return searchAlgorithm(pathMatrix, 0, 0, 3, 3);
    }



    public static void main(String[] args) {
        Scanner scan = new Scanner(System.in);
        String area = scan.nextLine();
        String it = scan.nextLine();
        System.out.println(findingMinimalPath(area, it));
    }
}
