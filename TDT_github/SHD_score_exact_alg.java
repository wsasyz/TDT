import java.io.*;
import java.util.*;

public class trial {
    static int[] from, to;
    static int moves;
    static int upperBound = 150;
    static double thre = 19.51142;// 18.18929
    static int depth = 5;
    static int[] contriB, contriC;
    static int[] dist;
    static int[][] counts;
    static void initSix() throws Exception{
        moves = depth * (depth + 1);
        from = new int[moves];
        to = new int[moves];
        int p = 0;
        for (int i = 0; i < depth + 1; i++)
            for (int j = 0; j < depth + 1; j++){
                if (i == j) continue;
                from[p] = i;
                to[p] = j;
                p++;
            }
        //(2,0),(0,2),(1,1),(1,0),(0,1). compute inner boundary
        contriB = new int[depth + 1];
        contriC = new int[depth + 1];
        contriB[0] = 2;
        contriB[2] = 1;
        contriB[3] = 1;
        contriC[1] = 2;
        contriC[2] = 1;
        contriC[4] = 1;
        counts = new int[depth + 1][upperBound + 1];
        for (int i = 0; i < depth + 1; i++)
            counts[i][0] = 1;
        for (int i = 0; i < upperBound + 1; i++)
            counts[0][i] = 1;
        for (int i = 1; i < depth + 1; i++)
            for (int j = 1; j < upperBound + 1; j++)
                counts[i][j] = counts[i - 1][j] + counts[i][j - 1];
        dist = new int[counts[depth][upperBound]];
    }
    static int findInd(int[] index){
        int uB = upperBound;
        int ind = 0;
        for (int i = 0; i < depth; i++){
            for (int j = 0; j < index[i]; j++){
                ind += counts[depth - i - 1][uB - j];
            }
            uB -= index[i];
        }
        return ind;
    }
    static int checkBoundary(int[] index){
        int b = 0, c = 0;
        for (int i = 0; i < depth; i++){
            b += contriB[i] * index[i];
            c += contriC[i] * index[i];
        }
        //if the difference between b and c is large, not inner boundary
        if (b + c != 0 && thre < (b-c)*(b-c)*1.0/(b+c))
            return 1;
        for (int i = 0; i < moves; i++){
            if (index[from[i]] == 0) continue;
            int b0, c0;
            b0 = b - contriB[from[i]] + contriB[to[i]];
            c0 = c - contriC[from[i]] + contriC[to[i]];
            if (b0 + c0 != 0 && thre < (b0 - c0) * (b0 - c0) * 1.0 / (b0 + c0))
                return 0;
        }
        return -1;
    }
    static void sixDim() throws Exception{
        int[] index, sums;
        int p, ind;
        Vector<int[]> vs = new Vector<int[]>(4000000);
        //iteration starts here
        index = new int[depth + 1];
        sums = new int[depth + 1];
        p = depth - 1;
        ind = 0;
        while (true){
            //insert code starts
            dist[ind] = (checkBoundary(index)) * 2 * upperBound;
            if (dist[ind] == 0) {
                vs.add(index.clone());
                if (vs.size() % 1000000 == 0){
                    System.out.println(vs.size());
                }
            }
            ind++;
            //inserted code ends
            index[p]++;
            while (p > 0 && sums[p] + index[p] == upperBound + 1) {
                sums[p]++;
                p--;
                index[p]++;
            }
            if (p == 0 && sums[0] + index[0] == upperBound + 1)
                break;
            for (;p != depth - 1; ){
                p++;
                index[p] = 0;
                sums[p+1] = sums[p];
            }
        }
        int round = 1;
        //iteration ends here
        Vector<int[]> vs2 = new Vector<int[]>();
        while (vs.size() != 0) {
            System.out.println("Round " + round + " pased!");
            round++;
            for (int[] vindex : vs){
                int curDist = dist[findInd(vindex)];
                for (int i = 0; i < moves; i++){
                    if (vindex[from[i]] == 0)
                        continue;
                    vindex[from[i]]--;
                    vindex[to[i]]++;
                    ind = findInd(vindex);
                    if (dist[ind] == 2 * upperBound){
                        dist[ind] = curDist + 1;
                        vs2.add(vindex.clone());
                    }
                    else if (dist[ind] == -2 * upperBound){
                        dist[ind] = curDist - 1;
                        vs2.add(vindex.clone());
                    }
                    vindex[from[i]]++;
                    vindex[to[i]]--;
                }
            }
            vs.clear();
            for (int[] vindex : vs2){
                vs.add(vindex);
            }
            vs2.clear();
        }
    }
    static void outputSix(BufferedWriter writer) throws Exception{
        StringBuilder sb = new StringBuilder();
        int[] index, sums;
        int p, ind;
        //iteration starts here
        index = new int[depth];
        sums = new int[depth + 1];
        p = depth - 1;
        ind = 0;
        while (true){
            //insert code starts
            /*for (int i = 0; i < depth; i++){
                sb.append(index[i]);
                sb.append(',');
            }*/
            sb.append(dist[ind]);
            sb.append('\n');
            if (sb.length() > 1000000) {
                writer.write(sb.toString());
                sb = new StringBuilder();
            }
            ind++;
            //inserted code ends
            index[p]++;
            while (p > 0 && sums[p] + index[p] == upperBound + 1) {
                sums[p]++;
                p--;
                index[p]++;
            }
            if (p == 0 && sums[0] + index[0] == upperBound + 1)
                break;
            for (;p != depth - 1; ){
                p++;
                index[p] = 0;
                sums[p+1] = sums[p];
            }
        }
        //iteration ends here
        writer.write(sb.toString());
    }
    static public void main2(String[] args) {
        try{
            initSix();
            BufferedReader reader = new BufferedReader(new FileReader("C:\\Users\\zhanglong\\Desktop\\tdt_input.txt"));
            Vector<String> inputs = new Vector<String>();
            Vector<Integer> inds = new Vector<Integer>(), inds2 = new Vector<Integer>();
            //reader.readLine();
            while(true){
                String s = reader.readLine();
                if (s == null) break;
                inputs.add(s);
                String[] line = s.split("[ \t,]");
                if (line.length != depth)
                    throw new Exception();
                int[] index = new int[depth];
                index[0] = Integer.parseInt(line[3]);
                index[1] = Integer.parseInt(line[4]);
                index[2] = Integer.parseInt(line[2]);
                index[3] = Integer.parseInt(line[0]);
                index[4] = Integer.parseInt(line[1]);
                int sum = index[0] + index[1] + index[2] + index[3] + index[4];
                int ind;
                if (sum > upperBound)
                    ind = -1;
                else
                    ind = findInd(index);
                inds.add(ind);
                inds2.add(ind);
            }
            reader.close();
            inds2.add(-1);
            Collections.sort(inds2);
            Hashtable<Integer, String> dict = new Hashtable<Integer, String>();
            reader = new BufferedReader(new FileReader("C:\\Users\\zhanglong\\Desktop\\dist" + upperBound + ".txt"));
            for (int i = 1; i < inds2.size(); i++){
                if (inds2.elementAt(i).equals(inds2.elementAt(i - 1))) continue;
                for (int j = 0; j < inds2.elementAt(i) - inds2.elementAt(i-1) - 1; j++){
                    reader.readLine();
                }
                String s = reader.readLine();
                if (s == null)
                    throw new Exception();
                //String[] line = s.split("[,\t ]");
                dict.put(inds2.elementAt(i), s);//line[line.length - 1]);
            }
            reader.close();
            BufferedWriter writer = new BufferedWriter(new FileWriter("output.txt"));
            StringBuilder sb = new StringBuilder();
            for (int i = 0; i < inputs.size(); i++){
                sb.append(inputs.elementAt(i));
                if (inds.elementAt(i) == -1)
                    sb.append(" NA\n");
                else{
                    sb.append(" ");
                    sb.append(dict.get(inds.elementAt(i)));
                    sb.append("\n");
                }
            }
            writer.write(sb.toString());
            writer.close();
        }
        catch (Exception e){
            e.printStackTrace();
        }
    }
    static public void main(String[] args){
        try{
            BufferedWriter writer = new BufferedWriter(new FileWriter("dist1.txt"));
            long start = System.currentTimeMillis();
            initSix();
            sixDim();
            System.out.println(System.currentTimeMillis() - start);
            outputSix(writer);
            writer.close();
        }
        catch (Exception e){
            e.printStackTrace();
        }
    }
}
