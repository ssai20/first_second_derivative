import java.util.*;
public class SortingNotOddsNumbers {

    public static int[] sortedNotOddsNumbers(int[] array)
    {
        if (array.length>10000) throw new RuntimeException("Quantity of numbers should be less than 10000.");
        List<Integer> arrayNotOddsNumbers = new ArrayList<>();
        for (int i=0;i<array.length; i++)
        {
            if (array[i]%2!=0)
            {
                arrayNotOddsNumbers.add(array[i]);
            }
        }
        Collections.sort(arrayNotOddsNumbers);
        int j=0;
        for (int i=0;i<array.length; i++)
        {
            if (array[i]%2!=0)
            {
                array[i] = arrayNotOddsNumbers.get(j);
                j++;
            }
        }
        return array;
    }

    public static void main(String[] args) {
        Scanner scan = new Scanner(System.in);
        String numbers = scan.nextLine();
        if (numbers.length()==0) System.out.println(numbers);
        if (numbers.length()>0)
        {
            int[] arrayNumbers = Arrays.stream(numbers.split(" ")).mapToInt(Integer::parseInt).toArray();
            StringBuilder stringBuilder = new StringBuilder();
            for (int i=0;i<arrayNumbers.length-1;i++)
            {
                stringBuilder.append(sortedNotOddsNumbers(arrayNumbers)[i]);
                stringBuilder.append(" ");
            }
            stringBuilder.append(sortedNotOddsNumbers(arrayNumbers)[arrayNumbers.length-1]);
            System.out.println(stringBuilder);
        }
    }
}
