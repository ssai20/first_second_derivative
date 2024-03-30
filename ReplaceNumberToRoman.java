public class ReplaceNumberToRoman{
    private static final String[] ROMAN_NUMERALS = {
            "M", "CM", "D", "CD", "C", "XC", "L", "XL", "X", "IX", "V", "IV", "I"
    };
    private static final int[] VALUES = {
            1000, 900, 500, 400, 100, 90, 50, 40, 10, 9, 5, 4, 1
    };

    public static String convertToRoman(int number) {
        if (number < 1 || number > 3999) {
            throw new IllegalArgumentException("Number out of range (1-3999)");
        }

        StringBuilder roman = new StringBuilder();

        for (int i = 0; i < VALUES.length; i++) {
            while (number >= VALUES[i]) {
                roman.append(ROMAN_NUMERALS[i]);
                number -= VALUES[i];
            }
        }

        return roman.toString();
    }

    public static void main(String[] args) {
        int number = 190;
        String romanNumeral = convertToRoman(number);
        System.out.println(romanNumeral);
    }
}
