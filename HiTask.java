import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;

public class HiTask {
    public static void main(String[] args) throws IOException {
        try (BufferedReader br = new BufferedReader(new InputStreamReader(System.in))) {
            int s = Integer.parseInt(br.readLine());
            if (s>7){
                System.out.println("Привет!");
            }
        }
    }
}
