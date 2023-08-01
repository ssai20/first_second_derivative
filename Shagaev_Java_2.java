//Составить алгоритм: если введенное имя совпадает с Вячеслав, то вывести “Привет, Вячеслав”,
// если нет, то вывести "Нет такого имени"

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;

public class Shagaev_Java_2 {
    public static void main(String[] args) {
        try (BufferedReader br = new BufferedReader(new InputStreamReader(System.in))) {
            String s = br.readLine();
            if (s.equals("Вячеслав")){
                System.out.println("Привет, Вячеслав");
            }
            else {
                System.out.println("Нет такого имени");
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
