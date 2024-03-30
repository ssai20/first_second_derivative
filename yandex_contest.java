import java.util.Scanner;

public class yandex_contest {
    public static void main(String[] args) {
        String s;
        Scanner scan = new Scanner(System.in);
        s = scan.nextLine();
        String[] stringArray = s.split(" ");

        System.out.println(Integer.parseInt(stringArray[0])+Integer.parseInt(stringArray[1]));
    }
}
