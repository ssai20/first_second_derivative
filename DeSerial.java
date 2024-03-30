import java.util.Base64;

public class DeSerial {
    public static void main(String[] args) {
        String encodedString = "SmF2YSDQvdC1INGC0L7RgNC80L7Qt9C40YIhCg==";
        byte[] decodedBytes = Base64.getDecoder().decode(encodedString);
        String decodedString = new String(Base64.getDecoder().decode(encodedString));
        System.out.println(decodedString);
    }
}