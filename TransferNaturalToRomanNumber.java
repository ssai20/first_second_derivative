import java.util.Scanner;

public class TransferNaturalToRomanNumber {
    public static String getOneToFour(int number)
    {
        StringBuilder oneToFour = new StringBuilder();
        if ( (number>=1)&&(number<=3) )
            for (int i=0;i<number;i++)
                oneToFour.append("I");
        if ( (number==4) )
            oneToFour.append("IV");
        return oneToFour.toString();
    }
    public static String getOneToNine(int number)
    {
        StringBuilder oneToNine = new StringBuilder();
        if ( (number>=1)&&(number<=4) )
            oneToNine.append(getOneToFour(number));
        if ( (number>=5)&&(number<=8) )
        {
            oneToNine.append("V");
            oneToNine.append(getOneToFour(number-5));
        }
        if ( (number==9) )
            oneToNine.append("IX");
        return oneToNine.toString();
    }
    public static String getOneToFortyNine(int number)
    {
        StringBuilder oneToFortyNine = new StringBuilder();
        if ( (number>=1) && (number<=9) )
            oneToFortyNine.append(getOneToNine(number));
        if ( (number>=10)&&(number<=39) )
        {
            for (int i = 0; i < (number - number % 10.) / 10.; i++)
                oneToFortyNine.append("X");
            oneToFortyNine.append(getOneToNine(number%10));
        }
        if ( (number>=40)&&(number<=49) )
        {
            oneToFortyNine.append("XL");
            oneToFortyNine.append(getOneToNine(number%10));
        }
        return oneToFortyNine.toString();
    }
    public static String getOneToNinetyNine(int number)
    {
        StringBuilder oneToNinetyNine = new StringBuilder();
        if ( (number>=1)&&(number<=49) )
            oneToNinetyNine.append(getOneToFortyNine(number));
        if ( (number>=50)&&(number<=89) )
        {
            oneToNinetyNine.append("L");
            oneToNinetyNine.append(getOneToFortyNine(number-50));
        }
        if ( (number>=90)&&(number<=99) )
        {
            oneToNinetyNine.append("XC");
            oneToNinetyNine.append(getOneToNine(number-90));
        }
        return oneToNinetyNine.toString();
    }
    public static String getOneToFourHundredNinetyNine(int number) {
        StringBuilder oneToFourHundredNinetyNine = new StringBuilder();
        if ( (number>=1)&&(number<100) )
            oneToFourHundredNinetyNine.append(getOneToNinetyNine(number));
        if ( (number>=100) && (number<=399))
        {
            for (int i=0;i<(number-number%100.)/100.;i++)
            {
                oneToFourHundredNinetyNine.append("C");
            }
            oneToFourHundredNinetyNine.append(getOneToNinetyNine(number%100));
        }
        if ((number>=400)&&(number<499))
        {
            oneToFourHundredNinetyNine.append("CD");
            oneToFourHundredNinetyNine.append(getOneToNinetyNine(number%100));
        }
        return oneToFourHundredNinetyNine.toString();
    }

    public static String getOneToNineHundredNineteenNine(int number)
    {
        StringBuilder oneToNineHundredNineteenNine = new StringBuilder();
        if ( (number>=1)&&(number<500) )
            oneToNineHundredNineteenNine.append(getOneToFourHundredNinetyNine(number));
        if ( (number>=500) && (number<=899) )
        {
            oneToNineHundredNineteenNine.append("D");
            oneToNineHundredNineteenNine.append(getOneToFourHundredNinetyNine(number-500));
        }
        if ( (number>=900) && (number<=999) )
        {
            oneToNineHundredNineteenNine.append("CM");
            oneToNineHundredNineteenNine.append(getOneToNinetyNine(number%100));
        }
        return oneToNineHundredNineteenNine.toString();
    }


    public static String getOneToThreeThousandNineHundredNineteenNine(int number)
    {
        StringBuilder oneToThreeThousandNineHundredNineteenNine = new StringBuilder();
        if ( (number>=1) && (number<=999) )
        {
            oneToThreeThousandNineHundredNineteenNine.append(getOneToNineHundredNineteenNine(number));
        }
        if ( (number>=1000) && (number<=3999) )
        {
            for (int i=0;i<(number-number%1000)/1000.;i++)
            {
                oneToThreeThousandNineHundredNineteenNine.append("M");
            }
            oneToThreeThousandNineHundredNineteenNine.append(getOneToNineHundredNineteenNine(number%1000));
        }

        return oneToThreeThousandNineHundredNineteenNine.toString();
    }


    public static void main(String[] args) {
        Scanner scan = new Scanner(System.in);
        int naturalNumber = scan.nextInt();
        System.out.println(getOneToThreeThousandNineHundredNineteenNine(naturalNumber));
    }
}
