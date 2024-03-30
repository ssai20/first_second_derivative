import java.util.*;

public class Main {
    public static void main(String[] args) {
        Scanner scanner = new Scanner(System.in);
        int n = scanner.nextInt();
        scanner.nextLine();
        if (n > 2.e5 || n < 1) {
            System.out.println("Необходимо ввести значения от 1 до 200000");
            return;
        }

        String inp = scanner.nextLine();
        String[] s = inp.split(" ");
        if (n != s.length) {
            System.out.println("Не соответствие количества чисел заданному.");
            return;
        }
        int[] a = new int[s.length];
        for (int i = 0; i < s.length; i++) {
            a[i] = Integer.parseInt(s[i]);
            if (a[i] > 2.e5 || a[i] < 1) {
                System.out.println("Необходимо ввести значения от 1 до 200000");
                return;
            }
        }

        HashMap<Integer, Integer> map = new HashMap<>();
        List<Integer> listL = new ArrayList<>();
        for (int i = 0; i < a.length; i++) {
            if (map.containsKey(a[i]))
                map.put(a[i], map.get(a[i]) + 1);
            else map.put(a[i], 1);


            List<Integer> list = new ArrayList<>(map.values());
            Collections.sort(list);
            list.get(list.size() - 1);

            if (((list.get(0) == list.get(list.size() - 1) - 1) && (list.get(0) == list.get(list.size() - 2)))
                    || ((list.get(0)+1 == list.get(list.size() - 1)) && (list.get(1) == list.get(list.size() - 1))))
            {

                listL.add(i+1);
            }
        }

        System.out.println(Collections.max(listL));


    }
}