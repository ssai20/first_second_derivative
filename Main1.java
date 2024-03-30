import java.util.Scanner;

public class Main1 {
    public static void main(String[] args) {
        Scanner scanner = new Scanner(System.in);

        String inp = scanner.nextLine();
        String[] heights = inp.split(" ");

        int h1 = Integer.parseInt(heights[0]);
        int h2 = Integer.parseInt(heights[1]);
        int h3 = Integer.parseInt(heights[2]);
        int h4 = Integer.parseInt(heights[3]);
        if (h1 > 300 || h1 < 0 || h2 > 300 || h2 < 0 || h3 > 300 || h3 < 0 ||h4 > 300 || h4 < 0 ) {
            System.out.println("Необходимо ввести значения от 0 до 300");
            return;
        }
        if ((h1<=h2&&h2<=h3&&h3<=h4) || (h1>=h2&&h2>=h3&&h3>=h4)) System.out.println("YES");
        else System.out.println("NO");
    }
}
