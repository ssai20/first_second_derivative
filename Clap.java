import java.util.Arrays;
import java.util.HashMap;
import java.util.Scanner;

public class Clap {

    public static int rowSearchNumber(int symbolNumbers, int[] textArray, HashMap mapIndexToRow)
    {
        int number = 0;
        int[] capacityArray = new int[symbolNumbers];
        for (int i=0;i<symbolNumbers;i++)
        {
            capacityArray[i] = (int)mapIndexToRow.get(textArray[i]);
        }
        for (int i=1;i<symbolNumbers;i++)
        {
            if (capacityArray[i]!=capacityArray[i-1])
            {
                number++;
            }
        }
        return number;
    }

    public static int textCapacity(int keysNumber, String keysIndex, String rowNumbers, int symbolNumbers, String text)
    {
        int capacity = 0;
        int[] keysIndexArray = Arrays.stream(keysIndex.split(" ")).mapToInt(Integer::parseInt).toArray();
        int[] rowNumbersArray = Arrays.stream(rowNumbers.split(" ")).mapToInt(Integer::parseInt).toArray();

        int[] textArray = Arrays.stream(text.split(" ")).mapToInt(Integer::parseInt).toArray();

        HashMap<Integer, Integer> mapIndexToRow = new HashMap<>();
        for (int i=0;i<keysNumber;i++)
        {
            mapIndexToRow.put(keysIndexArray[i], rowNumbersArray[i]);
        }
        capacity = rowSearchNumber(symbolNumbers, textArray, mapIndexToRow);
        return capacity;
    }


    public static void main(String[] args) {
        Scanner scan = new Scanner(System.in);
        int keysNumber = scan.nextInt();
        scan.nextLine();
        String keysIndex = scan.nextLine();
        String rowNumbers = scan.nextLine();
        int symbolNumbers = scan.nextInt();
        scan.nextLine();
        String text = scan.nextLine();
        System.out.println(textCapacity(keysNumber, keysIndex, rowNumbers, symbolNumbers, text));
    }
}
