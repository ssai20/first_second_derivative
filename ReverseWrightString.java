import java.util.Scanner;

public class ReverseWrightString {
    public static String reverseString(String simbolList)
    {
        StringBuilder reversedString = new StringBuilder();
        for (int i=simbolList.length()-1;i>=0;i--)
        {
            if (Character.isLetter(simbolList.charAt(i)))
                reversedString.append(simbolList.charAt(i));
        }
        return reversedString.toString();
    }
    public static void main(String[] args) {
        Scanner scan = new Scanner(System.in);
        String simbolList = scan.nextLine();
        System.out.println(reverseString(simbolList));
    }
}
//-d -i /Users/work/Desktop/out.txt /Users/work/Desktop/in1.txt /Users/work/Desktop/in2.txt /Users/work/Desktop/in1.txt /Users/work/Desktop/in3.txt