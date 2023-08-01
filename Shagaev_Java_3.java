import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.Arrays;

//Составить алгоритм: на входе есть числовой массив, необходимо вывести элементы массива кратные 3
public class Shagaev_Java_3 {
    public static void main(String[] args) {
        try (BufferedReader br = new BufferedReader(new InputStreamReader(System.in))) {
            String s = br.readLine();
            int[] array = Arrays.stream(s.split(" ")).mapToInt(Integer::parseInt).toArray();
            for (int elem:array){
                if (elem % 3 == 0)
                    System.out.print(elem + " ");
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
