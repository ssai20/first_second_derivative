import java.io.*;
import java.util.*;


class VariousTypes<T extends Comparable<T>> implements Comparable<T>{
    public T a;

    public int compareIntOrString(T a, T b)
    {
        if (a instanceof Integer) {
            System.out.println("Integer: ");
            return (Integer) a - (Integer) b;
        }
        if (a instanceof String){
            System.out.println("String: ");
            return a.compareTo(b);
        }
        return 0;
    }
    public List<T> genericMergeSortAlgorithm(List<T> listFile1, List<T> listFile2) {

//        if (a instanceof Integer) System.out.println("INTEGERRRR");
        List<T> listResult = new ArrayList<T>();
        if ((Objects.equals(listFile1, new ArrayList<T>()))&&(!Objects.equals(listFile2, new ArrayList<T>()))){
            return listFile2;
        }
        if ((!Objects.equals(listFile1, new ArrayList<T>()))&&(Objects.equals(listFile2, new ArrayList<T>()))){
            return listFile1;
        }
        if ((Objects.equals(listFile1, new ArrayList<T>()))&&(Objects.equals(listFile2, new ArrayList<T>()))){
            return listResult;
        }

        Queue<T> queueFile1 = new LinkedList<T>(listFile1);
        Queue<T> queueFile2 = new LinkedList<T>(listFile2);
        T buffer = null;
        boolean bufferInRightColumn = true;
//        if (queueFile1.peek()<queueFile2.peek()){
//        if (compareIntOrString(queueFile1.peek(), queueFile2.peek())<0) {
        if (queueFile1.peek().compareTo(queueFile2.peek()) < 0) {
//            if (queueFile1.peek() instanceof Integer) {
//                System.out.println("Integer: ");
//            }
            System.out.println(queueFile1.peek() + " < " + queueFile2.peek());
            buffer = queueFile2.peek();
            System.out.println("buffer = "+buffer);
            bufferInRightColumn = true;
            listResult.add(queueFile1.peek());
            queueFile1.poll();
//            queueFile2.poll();
        } else if (queueFile1.peek().compareTo(queueFile2.peek()) >= 0) {

            System.out.println(queueFile2.peek() + " < " + queueFile1.peek());
            bufferInRightColumn = false;
            buffer = queueFile1.peek();
            System.out.println("buffer = "+buffer);
            listResult.add(queueFile2.peek());
//            queueFile1.poll();
            queueFile2.poll();
        }

        while (!queueFile1.isEmpty() && !queueFile2.isEmpty()) {
            if (!queueFile1.isEmpty() || !queueFile2.isEmpty()) {
                if (bufferInRightColumn) {
//                    if (queueFile1.peek() instanceof Integer) {
//                        System.out.println("Integer: ");
//                    }
//                if (queueFile1.peek() < buffer) {
                    if (queueFile1.peek().compareTo(buffer) < 0) {
//                    if (compareIntOrString(queueFile1.peek(), buffer) < 0) {
                        listResult.add(queueFile1.peek());
                        System.out.println(queueFile1.peek() + " < " + buffer);
                        queueFile1.poll();
                        bufferInRightColumn = true;
                        System.out.println(bufferInRightColumn+"HHHH");
                    } else if (queueFile1.peek().compareTo(buffer) >= 0) {
                        System.out.println(buffer + " < " + queueFile1.peek());
                        listResult.add(buffer);
                        buffer = queueFile1.peek();
                        System.out.println("buffer = " + buffer);
                        queueFile2.poll();
                        bufferInRightColumn = false;
                        System.out.println(bufferInRightColumn);
                    }

                } else if (!bufferInRightColumn) {
//                if (queueFile2.peek() instanceof Integer) {
//                    System.out.println("Integer: ");
//                }
//                if (queueFile2.peek() < buffer) {
                if (queueFile2.peek().compareTo(buffer) < 0) {
//                    if (compareIntOrString(queueFile2.peek(), buffer) < 0) {
                        listResult.add(queueFile2.peek());
                        System.out.println(queueFile2.peek() + " < " + buffer);
                        System.out.println("buffer = " + buffer);
                        queueFile2.poll();
                        bufferInRightColumn = false;
                        System.out.println(bufferInRightColumn);
                    } else if (queueFile2.peek().compareTo(buffer) >= 0){
                        System.out.println(buffer + " < " + queueFile2.peek());

                        listResult.add(buffer);
                        buffer = queueFile2.peek();
                        System.out.println("buffer = " + buffer);
                        queueFile1.poll();
                        bufferInRightColumn = true;
                        System.out.println(bufferInRightColumn);
                    }
                }
            }


            if (queueFile1.isEmpty() && !queueFile2.isEmpty()) {
                System.out.println("Левая колонка пуста");
                listResult.add(buffer);
                queueFile2.poll();
                listResult.addAll(queueFile2);
                queueFile2.removeAll(queueFile2);
            }
            if (!queueFile1.isEmpty() && queueFile2.isEmpty()) {
                System.out.println("Правая колонка пуста");
                listResult.add(buffer);
                queueFile1.poll();
                listResult.addAll(queueFile1);
                queueFile1.removeAll(queueFile1);
            }
        }

        for(int k = 0;k<listResult.size();k++)
        {
            System.out.println(listResult.get(k));
        }

        return listResult;

}

    public List<T> mapFileToList(File file, String dataType)
    {
        List<T> listFile = new ArrayList<>();

        try (BufferedReader br = new BufferedReader(new InputStreamReader(new FileInputStream(file), "UTF-8"))){

            while (true)
            {
                String string = br.readLine();
                if (string == null) {
                    break;
                }
                if (dataType.equals("-i")) {
                    Integer i = Integer.parseInt(string);
                    listFile.add((T) i);
                }
                if (dataType.equals("-s")) {
                    String i = string;
                    listFile.add((T) i);
                }
            }
        } catch (IOException e1){
            e1.printStackTrace();
        }
        return listFile;
    }

    public void mapListToFile (String fileOutPath, List<T> resultList)
    {
        File file = new File(fileOutPath);
        try (BufferedWriter bw = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(file),  "UTF-8"))){
            for (int i=0;i<resultList.size();i++) {
                bw.write(resultList.get(i).toString());
                bw.newLine();
            }
        } catch (IOException e1){
            e1.printStackTrace();
        }
    }

    @Override
    public int compareTo(T o) {
        if ((this.a instanceof Integer) && (o instanceof Integer)){

            return (Integer)this.a-(Integer)o;
        }
        if ((this.a instanceof String) && (o instanceof String)){
            return ((String) this.a).compareTo((String) o);
        }
        return 0;
    }

//    private int compare(T a, T b)
//    {
//        return a.compareTo(b);
//    }
//    @Override
//    public int compareTo(T o) {
//        int result;
////        result = Double.compare(d, o.d);
////        if(result != 0) return result;
//        result = Integer.compare((Integer) a, (Integer) o);
//        if(result != 0) return result;
////        result = s.compareTo(o.s);
////        if(result != 0) return result;
////        result = Character.compare(c, o.c);
//        return result;
//    }
//    @Override

//    public int compareTo(T o) {
////        if ((this.a instanceof Integer) && (o instanceof Integer)){
////            System.out.println("INTEGER");
////            String w = a.toString();
////            String q = o.toString();
////            int a1 = Integer.parseInt(w);
////            int a2 = Integer.parseInt(q);
////            if (a1<a2) {
////                System.out.println(a1+"<<"+a2);
////                return -1;
////            }
////            if (a1>a2) return 1;
////            if (a1==a2) return 0;
//            return (Integer)this.a-(Integer)o;
//            //return ((Integer) this.a).compareTo((Integer) o);
////        }
////        if ((this.a instanceof String) && (o instanceof String)){
////            return ((String) this.a).compareTo((String) o);
////        }
////        return 0;
//    }

}





public class MergeSort{
//
//
//    public static List<Integer> mergeSortAlgorithm(List<Integer> listFile1, List<Integer> listFile2){
//        List<Integer> listResult = new ArrayList<>();
//        Queue<Integer> queueFile1 = new LinkedList<>(listFile1);
//        Queue<Integer> queueFile2 = new LinkedList<>(listFile2);
//        int buffer;
//        boolean bufferInRightColumn;
//        if (queueFile1.peek()<queueFile2.peek()){
//            buffer = queueFile2.peek();
//            bufferInRightColumn = true;
//            listResult.add(queueFile1.peek());
//            queueFile1.poll();
//            queueFile2.poll();
//        }
//        else {
//            bufferInRightColumn = false;
//            buffer = queueFile1.peek();
//            listResult.add(queueFile2.peek());
//            queueFile1.poll();
//            queueFile2.poll();
//        }
//
//        while (!queueFile1.isEmpty()||!queueFile2.isEmpty())
//        {
//            if (!queueFile1.isEmpty()&&!queueFile2.isEmpty()) {
//                if (bufferInRightColumn) {
//                    if (queueFile1.peek() < buffer) {
//                        listResult.add(queueFile1.peek());
//                        System.out.println(queueFile1.peek() + " < " + buffer);
//                        queueFile1.poll();
//                        bufferInRightColumn = true;
//                    } else {
//                        System.out.println(buffer + " < " + queueFile1.peek());
//                        listResult.add(buffer);
//                        buffer = queueFile1.peek();
//                        queueFile1.poll();
//                        bufferInRightColumn = false;
//                    }
//
//                } else {
//                    if (queueFile2.peek() < buffer) {
//                        listResult.add(queueFile2.peek());
//                        System.out.println(queueFile2.peek() + " < " + buffer);
//                        queueFile2.poll();
//                        bufferInRightColumn = false;
//                    } else {
//                        System.out.println(buffer + " < " + queueFile2.peek());
//                        listResult.add(buffer);
//                        buffer = queueFile2.peek();
//                        queueFile2.poll();
//                        bufferInRightColumn = true;
//                    }
//                }
//            }
//
//            if (queueFile1.isEmpty()&&!queueFile2.isEmpty()) {
//                listResult.add(buffer);
//                listResult.addAll(queueFile2);
//                queueFile2.removeAll(queueFile2);
//            }
//            if (!queueFile1.isEmpty()&&queueFile2.isEmpty()) {
//                listResult.add(buffer);
//                listResult.addAll(queueFile1);
//                queueFile1.removeAll(queueFile1);
//            }
//        }
//        for(int k=0;k<listResult.size();k++)
//        {
//            System.out.println(listResult.get(k));
//        }
//        return listResult;
//    }





    public static void mergeSort(String conditions) {

        String[] s = conditions.split(" ");
        if (!s[0].equals("-a") && !s[0].equals("-d")) {
            StringBuilder conditionsSB = new StringBuilder(conditions);
            conditionsSB.insert(0, "-a ");
            conditions = conditionsSB.toString();
            System.out.println(conditions);
        }
        List<String> conditionList = List.of(conditions.split(" "));
        String sortType = conditionList.get(0);
        String dataType = conditionList.get(1);
        String fileOutPath = conditionList.get(2);

        if (dataType.equals("-i")) {
            VariousTypes<Integer> variousTypes = new VariousTypes<>();
            List<Integer> resultList = new ArrayList<>();


            for (int i = 3; i < conditionList.size(); i++) {
                resultList = variousTypes.genericMergeSortAlgorithm(resultList, variousTypes.mapFileToList(new File(conditionList.get(i)), dataType));
            }
            if (sortType.equals("-a")) {
                variousTypes.mapListToFile(fileOutPath, resultList);
            }
            if (sortType.equals("-d")) {
                Collections.reverse(resultList);
                variousTypes.mapListToFile(fileOutPath, resultList);
            }
    }

        if (dataType.equals("-s")) {
            VariousTypes<String> variousTypes = new VariousTypes<>();
            List<String> resultList = new ArrayList<>();
            for (int i = 3; i < conditionList.size(); i++) {
                resultList = variousTypes.genericMergeSortAlgorithm(resultList, variousTypes.mapFileToList(new File(conditionList.get(i)), dataType));
            }
            if (sortType.equals("-a")) {
                variousTypes.mapListToFile(fileOutPath, resultList);
            }
            if (sortType.equals("-d")) {
                Collections.reverse(resultList);
                variousTypes.mapListToFile(fileOutPath, resultList);
            }
        }


//        VariousTypes<Integer> variousTypes = new VariousTypes<>();
//        List<Integer> resultList = sortAlgorithm.genericMergeSortAlgorithm(listFile1, listFile2);

//        List<Integer> listFile1 = new ArrayList<>();
//
//        try (BufferedReader br = new BufferedReader(new InputStreamReader(new FileInputStream(file1), "UTF-8"))){
//
//            while (true)
//            {
//                String string = br.readLine();
//                if (string == null) {
//                    break;
//                }
//
//                int i = Integer.parseInt(string);
//                listFile1.add(i);
//            }
//        } catch (IOException e1){
//            e1.printStackTrace();
//        }
//
//        List<Integer> listFile2 = new ArrayList<>();
//        try (BufferedReader br = new BufferedReader(new InputStreamReader(new FileInputStream(file2), "UTF-8"))){
//
//            while (true)
//            {
//                String string = br.readLine();
//                if (string == null) {
//                    break;
//                }
//                int j = Integer.parseInt(string);
//                listFile2.add(j);
//            }
//        } catch (IOException e1){
//            e1.printStackTrace();
//        }


//        VariousTypes<Integer> sortAlgorithm = new VariousTypes<>();
//        List<Integer> resultList = sortAlgorithm.genericMergeSortAlgorithm(listFile1, listFile2);
//        for (int i=0;i< inFile.length;i++){
//
//        }

//        try (BufferedWriter bw = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(new File(fileOutPath)),  "UTF-8"))){
//            for (int i=0;i<resultList.size();i++) {
//                bw.write(resultList.get(i).toString());
//                bw.newLine();
//            }
//        } catch (IOException e1){
//            e1.printStackTrace();
//        }


    }
    public static void main(String[] args) {
        String conditions = null;
        try (BufferedReader br = new BufferedReader(new InputStreamReader(System.in))) {
//            System.out.println("Путь к первому файлу: ");
//            filePath1 = br.readLine();
//            System.out.println("Путь ко второму файлу: ");
//            filePath2 = br.readLine();
            System.out.println("Введите условия сортировки и входные и входной файл: ");
            conditions = br.readLine();
        } catch (IOException e) {
            e.printStackTrace();
        }

        mergeSort(conditions);
//        mergeSort(filePath1, filePath2, resultFilePath);



//        mergeSort("/Users/work/Desktop/in1.txt", "/Users/work/Desktop/in2.txt", "/Users/work/Desktop/out.txt");
    }
}
//абсолютный путь к файлам
//-d -i /Users/work/Desktop/out.txt /Users/work/Desktop/in1.txt /Users/work/Desktop/in2.txt

//  -i ./src/main/java/out.txt ./src/main/java/in1.txt ./src/main/java/in2.txt

//class MessageReaderFromFile{
//
//    public static List<LocalMessage> readAllMessagesFromTextFile() throws IOException
//    {
//        String fileName = "/Users/work/Desktop/emails.txt";
//        File file = new File(fileName);
//        List<LocalMessage> messagesSync = Collections.synchronizedList(new ArrayList<>());
//        try (BufferedReader br = new BufferedReader(new InputStreamReader(new FileInputStream(file), "UTF-8"))){
//            while (true)
//            {
//                String email = br.readLine();
//                if (email!=null) {
//                    LocalMessage message = new LocalMessage(email);
//                    messagesSync.add(message);
//                }
//                if (email==null) {
//                    return messagesSync;
//                }
//
//            }
//        } catch (UnsupportedEncodingException e1){
//            e1.printStackTrace();
//        }
//        return messagesSync;
//    }
//}