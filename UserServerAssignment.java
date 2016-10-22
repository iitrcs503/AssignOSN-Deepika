package userserverassignment;
import java.io.*;
import java.math.RoundingMode;
import java.text.DecimalFormat;
import java.util.*;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author Deepika,Monika,Sneha
 */
public class UserServerAssignment {

    //declaring variables
    static int n ;   //No of users
    static int m = 0;   //No of servers
    static int k = 0;   //No of max replicas
    static int nmax;    //max users supported by system
    static int thresh = 45;

    static List<Integer> rel[]; //Array of list of users who have relationship among each other

    static double cost_r = 0.01;    //Average cost of reading
    static double cost_w = 1;    //Average cost of writing
    static double cost_t = 0.02;    //Average cost of transferring

    static int P[][];   //Assignment matrix
    static int X[][];   //Replication matrix

    static double avg_r = 1;    //Avgerage rate of reading
    static double avg_w = 0.1;  //Average rate of writing

    static int K[];    //No of replicas of each user

    static double load[];   //Load at each server
    
    static char mode = 'o'; //mode of algo
    
    static double imbalance[];
    
    static double alpha = 0.5;
    static double beta = 0.5;
    static double prob_reassign;
    static double prob_replicate;
    static double prob_remove;
    static double prob_threshold = 1.0;
    static double Tb;
    
    static Server server[];
    
    static BufferedWriter bw;
    
    public static void main(String[] args) {

       //take user input for m,k 
       Scanner sc = new Scanner(System.in);
       System.out.println("Enter number of users:");
       n = sc.nextInt();
       System.out.println("Enter number of servers:");
       m = sc.nextInt(); 
       System.out.println("Enter number of replicas:");
       k = sc.nextInt();
       nmax = m * thresh;
       P = new int[nmax][m];
       X = new int[nmax][m];
       K = new int[nmax];
       load = new double[m];
       rel = new ArrayList[nmax];
       imbalance = new double[m];

       int i,j;
       
        try {
            //take file input and initialize rel list
            BufferedReader br = new BufferedReader(new FileReader("input.txt"));
            bw = new BufferedWriter(new FileWriter("output.txt"));
            for(i=0;i<nmax;i++){
                rel[i] = new ArrayList<>();
            }
            String line;
            while ((line = br.readLine()) != null) {
               String x[] = line.split(" ");
               int u1 = Integer.parseInt(x[0]);
               int u2 = Integer.parseInt(x[1]);
               rel[u1].add(u2);
               rel[u2].add(u1);
            }

            //initialize replication matrix
            for(i=0;i<n;i++){
                for(j=0;j<m;j++){
                    X[i][j] = 0;
                }
            }   
            //initialize assignment matrix
            for(i =0; i < m; i++){
                for(j = 0; j < n; j++){
                    int value = (j%m);
                    if(value == i){
                        P[j][i] = 1;
                    }
                }
            }

            //create instances for server and initialize n, nr, V, Vr for each server
            server = new Server[m];
            int countUser;
            for(int a = 0; a < m; a++){
                countUser = 0;
                server[a] = new Server();
                server[a].V = new ArrayList<>();
                server[a].Vr = new ArrayList<>();
                for(int b = 0; b < n; b++){
                    if(P[b][a] == 1){
                        countUser++;
                        server[a].V.add(b);
                    }
                }
                server[a].n = countUser;
                server[a].nr = 0;
            }

            //calculate load of every server
            for(i=0;i<m;i++){
                server[i].L = readTransferLoad(server[i].n, server[i].V, i);
                double r_load = readLoad(server[i].n, server[i].V, i, server[i].L);
                double w_load = writeLoad(server[i].n, i, server[i].L, server[i].nr);
                load[i] = r_load + w_load;   
            }

            printResult();

            mainAlgo();
            addUser();
            br.close();
            bw.close();
        } catch (FileNotFoundException ex) {
            Logger.getLogger(UserServerAssignment.class.getName()).log(Level.SEVERE, null, ex);
        } catch (IOException ex) {
            Logger.getLogger(UserServerAssignment.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
    
    //ALGORITHM 1 
    public static void mainAlgo() throws IOException{
        mode = 'o';
        int x = 0;
        do{ 
            if(mode == 'o'){ 
                bw.write("**********Optimization Mode***********\n");
                System.out.println("**********Optimization Mode**********");
                boolean action = Optimization();
                printResult();
                if(action == false){
                    mode = 'c';
                }
            }
            else{
                bw.write("***********Convergence Mode***********\n");
                System.out.println("***********Convergence Mode***********");
                Convergence();
                printResult();
                mode = 'o';
                x=1;
            }
        }while(x == 0);
    }
    
    public static double[][] readTransferLoad(int nu,List<Integer> V,int id) {
        double l[][]=new double[nmax][m];
        
        for(int p = 0;p<n;p++){
            for(int g = 0;g<m;g++){
                l[p][g] = -1;
            }
        }
        for(int i=0;i<V.size();i++) {// n is users on a server
            double sum = 0;
            int uid = 0;
            int user = V.get(i);   //V is list of users on a server ,user is userid
            int len = (rel[user]).size();  //len consists of no of users related to user
            for(int t=0;t<m;t++){
                if(t == id){
                    continue;
                }
                        
                for(int j=0;j<len;j++){
                    uid = rel[user].get(j);
                    sum += P[uid][t] * (1 - X[user][t]) * cost_r * avg_r;    
                }
                l[user][t] = sum;
            }
        }
        return l;     
    }

    //function to calculate read load
    public static double readLoad(int nu,List<Integer> V,int id,double l[][]){
        double read_load = 0;
        for(int u=0; u<V.size();u++){  
           double sum = 0;
           int user = V.get(u);
            for(int i=0;i<m;i++){
                if(i == id)
                    continue;
                if(l[user][i] != -1)        
                    sum += l[user][i];
            }
            read_load += (sum + (avg_r * cost_r));
        }
        return read_load;  
    }
    
    //function to calculate write load
    public static double writeLoad(int n,int id,double l[][],int nr){
        double w1_load = 0,w2_load = 0,write_load = 0;
        for(int u=0;u<n;u++){  
            w1_load += (avg_w * cost_w);
        }
        for(int u=0;u<nr;u++){  
            w2_load += (avg_w * cost_w);
        }
        write_load = w1_load + w2_load;
        return write_load;
    }
    
    public static boolean Optimization() throws IOException {
        boolean action = false;
        int o = overloadServer(load);
        
        for(int i=0;i<server[o].V.size();i++){
            
            int user = server[o].V.get(i);
            int minLoaded = -1;
            double F_Reassign = Double.POSITIVE_INFINITY; 
            for(int s=0;s<m;s++){
                if(o != s){
                    double val = FReassign(server[o].L, o, s, user, load);
                    if(val < F_Reassign){
                        F_Reassign = val;
                        minLoaded = s;
                    }        
                }
            }

            double F_Replicate = Double.POSITIVE_INFINITY;
            int maxRel = -1;
            if(K[user]<k){
                for(int s=0;s<m;s++){
                    if(o != s){
                        double val = FReplica(server[o].L, user, s, load);
                        if(val < F_Replicate){
                            F_Replicate = val;
                            maxRel = s;
                        }        
                    }
                }
            }else{
                F_Replicate = Double.POSITIVE_INFINITY;
            }
            if(F_Reassign >= 0 && F_Replicate >= 0){
            }else if(F_Reassign < 0 && F_Replicate > 0){
                changeDueToReassign(user, o, minLoaded);
                i = -1;
                action = true;
            }else if(F_Reassign > 0 && F_Replicate < 0){
                changeDueToReplicate(user, maxRel, o);
                action = true;
            }else if(F_Reassign <= F_Replicate){
                changeDueToReassign(user, o, minLoaded);
                i = -1;
                action = true;
            }else if(F_Reassign > F_Replicate){
                changeDueToReplicate(user, maxRel, o);
                action = true;
            }
              
        }
        
        for(int i=0;i<server[o].Vr.size();i++){
            
            int user = server[o].Vr.get(i);
            //find primary server of user whose replica is removed from overloaded server
            int priServer = findUserPrimaryServer(user);
            double F_Remove = FRemoving(user, o, load);
            if(F_Remove < 0){
                changeDueToReplicaRemoval(priServer, o, user);
                action = true;
            }
        }    
        return action;
    }
    
    //function to find overloaded server
    public static int overloadServer(double load[]) throws IOException{
        double max = load[0];
        int idx = 0;
        for(int i=1; i<m; i++){
            if(max < load[i]){
                max = load[i];
                idx = i;
            }
        }
        bw.write("Overloaded server : "+idx+"\n");
        System.out.println("Overloaded server : "+idx);
        return idx;    
    }
    
    public static int minLoadedServer(double load[], int user) throws IOException{
        double min = 100000;
        int idx = 0;
        for(int i=0; i<m; i++){
            if(P[user][i] != 1){
                if(min > load[i]){
                    min = load[i];
                    idx = i;
                }
            }
            
        }
        bw.write("Minimum loaded server : "+idx+"\n");
        System.out.println("Minimum loaded server : "+idx);
        return idx;    
    }
    
    //Function to find server to which replicate user
    public static int userReplica(double L[][],int u) throws IOException{
        double max = -2;
        int s = 0;
        for(int i=0;i<m;i++){
            if(X[u][i] != 1){
                if(max < L[u][i]){
                    max = L[u][i];
                    s = i;
                }
            }            
        }
        bw.write("Server on which to replicate : "+s+"\n");
        System.out.println("Server on which to replicate : "+s);
        return s;
    }
    
    //function to find server to remove replica of user
    public static int findUserPrimaryServer(int u){ 
        for(int a=0;a<m;a++){
            if(P[u][a]==1)
                break;
        }
        return k;  
    }
    
    public static double FReassign(double L[][],int o, int minLoaded, int v,double load[]) throws IOException{
        double Fre;
        if(server[minLoaded].n < thresh){
            double delta_load = 0;
            for(int i=0;i<rel[v].size();i++){
                int uid = rel[v].get(i);
                delta_load += (P[uid][o] * avg_r * cost_r ) ;//
            }
            delta_load = delta_load - L[v][minLoaded] - ( X[v][minLoaded] * avg_w * cost_w );  
            double total_load = totalLoad(load);
            double new_load = total_load + delta_load;
            double F = calculateF(total_load);
            double F_new = calculateF(new_load);
            Fre = F_new - F;
        }
        else{
            Fre = Double.POSITIVE_INFINITY;
        }
        bw.write("Reassign user "+v+" from server "+o+" to server "+minLoaded+"? FReassign value : "+Fre+"\n");
        System.out.println("Reassign user "+v+" from server "+o+" to server "+minLoaded+"? FReassign value : "+Fre);
        return Fre;
         
    }
    
    public static double FReplica(double L[][],int v,int id, double load[]) throws IOException{
        double dp_load = avg_w * cost_w - L[v][id];
        double tp_load = totalLoad(load);
        double np_load = tp_load + dp_load;
        double Fp = calculateF(tp_load);
        double F_newp = calculateF(np_load);
        double Frep = F_newp - Fp;
        bw.write("Create replica of user "+v+" on server "+id+"? FReplica value : "+Frep+"\n");
        System.out.println("Create replica of user "+v+" on server "+id+"? FReplica value : "+Frep);
        return Frep;
    }
           
    public static double FRemoving(int v,int id, double load[]) throws IOException{
        double db_load = 0;
        for(int i=0;i<rel[v].size();i++){
            int uid = rel[v].get(i);
            db_load += (P[uid][id] * avg_r * cost_r);
        }
        db_load -= (cost_w * avg_w);
        double tb_load = totalLoad(load);
        double nb_load = tb_load + db_load;
        double F = calculateF(tb_load);
        double F_newb = calculateF(nb_load);
        double Fb = F_newb - F;
        bw.write("Remove replica of user "+v+" from server "+id+"? Fremove value : "+Fb+"\n");
        System.out.println("Remove replica of user "+v+" from server "+id+"? Fremove value : "+Fb);
        return Fb;

    }
    
    //function to change server due to reassign
    public static void changeDueToReassign(int user,int o,int minLoaded) throws IOException{
        bw.write("###### Reassigning user "+user+" from server "+o+" to server "+minLoaded+" ######\n");
        System.out.println("###### Reassigning user "+user+" from server "+o+" to server "+minLoaded+" ######");
        server[o].n -= 1;
        server[minLoaded].n += 1;
        server[o].V.remove(new Integer(user));
        server[minLoaded].V.add(user);
        P[user][o] = 0;
        P[user][minLoaded] = 1;
        for(int i=0;i<m;i++){
            server[o].L[user][i] = -1;
        }
        server[minLoaded].L = readTransferLoad(server[minLoaded].V.size(),server[minLoaded].V,minLoaded);
        double r_load = readLoad(server[o].n, server[o].V, o, server[o].L);
        double w_load = writeLoad(server[o].n, o, server[o].L, server[o].nr);
        load[o] = r_load + w_load;
        r_load = readLoad(server[minLoaded].n, server[minLoaded].V,minLoaded , server[minLoaded].L);
        w_load = writeLoad(server[minLoaded].n, minLoaded, server[minLoaded].L, server[minLoaded].nr);
        load[minLoaded] = r_load + w_load; 
        bw.write("Changes due to REASSIGNMENT :\n");
        System.out.println("Changes due to REASSIGNMENT :");
        printServer(o);
        printServer(minLoaded);
        
    }
    
    //function to change server due to replication
    public static void changeDueToReplicate(int user, int maxrel,int o) throws IOException{
        bw.write("###### Creating replica of user "+user+" on server "+maxrel+" ######\n");
        System.out.println("###### Creating replica of user "+user+" on server "+maxrel+" ######");
        server[maxrel].Vr.add(user);
        server[maxrel].nr += 1;
        K[user] += 1;
        X[user][maxrel] = 1;
        server[o].L = readTransferLoad(server[o].V.size(),server[o].V,o);
        double r_load = readLoad(server[o].n, server[o].V, o, server[o].L);
        double w_load = writeLoad(server[o].n, o, server[o].L, server[o].nr);
        load[o] = r_load + w_load;
        r_load = readLoad(server[maxrel].n, server[maxrel].V,maxrel , server[maxrel].L);
        w_load = writeLoad(server[maxrel].n, maxrel, server[maxrel].L, server[maxrel].nr);
        load[maxrel] = r_load + w_load; 
        bw.write("Changes due to REPLICATION :\n");
        System.out.println("Changes due to REPLICATION :");
        printServer(o);
        printServer(maxrel);
        
    }
    
    //function to change server dur to replica removal
    public static void changeDueToReplicaRemoval(int priServer,int o,int user) throws IOException{
        bw.write("###### Removing replica of user "+user+" from server "+o+" ######\n");
        System.out.println("###### Removing replica of user "+user+" from server "+o+" ######");
        server[o].nr -= 1;
        server[o].Vr.remove(new Integer(user));
        X[user][o] = 0;
        K[user] -= 1;
        server[priServer].L = readTransferLoad(server[priServer].V.size(),server[priServer].V,priServer);
        double r_load = readLoad(server[o].n, server[o].V, o, server[o].L);
        double w_load = writeLoad(server[o].n, o, server[o].L, server[o].nr);
        load[o] = r_load + w_load;
        r_load = readLoad(server[priServer].n, server[priServer].V,priServer , server[priServer].L);
        w_load = writeLoad(server[priServer].n, priServer, server[priServer].L, server[priServer].nr);
        load[priServer] = r_load + w_load; 
        bw.write("Changes due to REPLICA REMOVAL :\n");
        System.out.println("Changes due to REPLICA REMOVAL :");
        printServer(o);
        printServer(priServer);
        
    }

    //function to calculate total load
    public static double totalLoad(double load[]) throws IOException{
        double total_load = 0; 
        for(int i = 0; i<m;i++){   
            total_load += load[i];   
        }
        return total_load;    
    }
    
    //function to calculate F
    public static double calculateF(double total_load){
        
        double avg_load = total_load/m;
        for(int i = 0; i<m;i++){   
            imbalance[i] = avg_load - load[i];
        }

        double F = total_load + alpha * findMax(imbalance);
        return F;
    }
    
    //function to get maximum imbalance load
    public static double findMax(double im[])
    {
        double max = im[0];
        for(int j=1; j<m ;j++){
            if(max < im[j])
                max = im[j];
        }
        return max;
    }
     
    //Algorithm with convergence mode
    public static void Convergence() throws IOException
    {
        //Tb is a function of global objective value
        Tb = calculateTb();
        for(int i = 0; i < m; i++){
            for(int j = 0; j < server[i].V.size(); j++){
                int v = server[i].V.get(j);
                boolean reassigned = false;
                for(int k1 = 0; k1 < m; k1++){
                    if(i != k1) {
                        //Compute Reassign
                        if(!reassigned){
                            double F_Reassign = FReassign(server[i].L, i, k1, v, load);
                            if(F_Reassign == 0)
                                prob_reassign = 0;
                            else
                                prob_reassign = Math.min(1, java.lang.Math.pow(Math.E, (-F_Reassign/Tb)));
                            bw.write("Probability of reassigning user "+v+" from server "+i+" to server "+k1+" : "+prob_reassign+"\n");
                            System.out.println("Probability of reassigning user "+v+" from server "+i+" to server "+k1+" : "+prob_reassign);
                            if(prob_reassign >= prob_threshold){
                                reassigned = true;
                                changeDueToReassign(v, i, k1);
                            }
                        }
                        if(X[v][k1] == 0 && P[v][k1] != 1 && K[v] < k){
                            //Compute FRep
                            double F_Replicate = FReplica(server[i].L, v, k1, load);
                            if(F_Replicate == 0)
                                prob_replicate = 0;
                            else
                                prob_replicate = Math.min(1, java.lang.Math.pow(Math.E, (-F_Replicate/Tb)));
                            bw.write("Probability of replicating user "+v+" on server "+k1+" : "+prob_replicate+"\n");
                            System.out.println("Probability of replicating user "+v+" on server "+k1+" : "+prob_replicate);
                            if(prob_replicate >= prob_threshold){
                                changeDueToReplicate(v, k1, i);
                            }
                        }
                    }
                }
            }
            for(int l = 0; l < server[i].Vr.size(); l++){
                int user = server[i].Vr.get(l);
                //Compute F remove
                double F_Remove = FRemoving(user, i, load);
                if(F_Remove == 0)
                    prob_remove = 0;
                else
                    prob_remove = Math.min(1, java.lang.Math.pow(Math.E, (-F_Remove/Tb)));
                bw.write("Probability of removing replica of user "+user+" from server "+i+" : "+prob_remove+"\n");
                System.out.println("Probability of removing replica of user "+user+" from server "+i+" : "+prob_remove);
                if(prob_remove >= prob_threshold){
                    int priServer = findUserPrimaryServer(user);
                    changeDueToReplicaRemoval(priServer, i, user);
                }
            }    
        }
    }
    
    //function to calculate Tb
    public static double calculateTb() throws IOException{
        double F = 0, FL = 0, FU = 0;
        double total_load = totalLoad(load);
        F = calculateF(total_load);
        FL = calculateFL();
        FU = calculateFU(FL);
        Tb = beta *((F - FL)/(FU - FL));
        System.out.println("Value of Tb : "+Math.abs(Tb));
        return Math.abs(Tb);
    }
    
    //function to calculate FL
    public static double calculateFL(){
        double w1_load = 0,w2_load = 0, FL = 0;
        for(int i = 0; i < n; i++){
            w1_load = (avg_w * cost_w);
            w2_load = (avg_r * cost_r);
            int len = rel[i].size();  //len consists of no of users related to user
            for(int j = 0; j < len; j++){
                FL += (avg_r * cost_r);
            }
            
            FL += w1_load + w2_load; 
        }
        return FL;
    }
    
    //function to calculate FU
    public static double calculateFU(double FL){
        double FU;
        FU = FL + Math.max(calculateLc_max(), calculateLrep_max()) + ((alpha *(n - 1)* FL)/n);
        return FU;
    }
    
    //function to calculate Lc_max for FU
    public static double calculateLc_max() {
        double Lc_max = 0;
        for(int i = 0; i < n; i++){
            int len = rel[i].size();
            for(int j = 0; j < len; j++){
                for(int k = 0; k < m; k++){
                    Lc_max += (avg_r * cost_t)*(1 - P[rel[i].get(j)][k]);
                }
            }
        }
        return Lc_max;
    }
    
    //function to calculate Lrep_max for FU
    public static double calculateLrep_max(){
        double Lrep_max = 0;
        for(int i = 0; i < n; i++){
            int len = rel[i].size();
            Lrep_max += (avg_w * cost_w) * len;
        }
        return Lrep_max;
    }   
    
    //function to add user
    public static void addUser() throws IOException{
        int ch,ch1;
        do{
            if(n==nmax){
                bw.write("Max no. of users reached\n");
                System.out.println("Max no. of users reached");
                break;
            }
            System.out.println("Want to add user : 0/1");
            Scanner sc = new Scanner(System.in);
            ch = sc.nextInt();
            if(ch==1){
                n+=1;
                ch1 = n-1;
                bw.write("Adding user "+ch1+"\n");
                System.out.println("Adding user "+ch1);
                createRelation(ch1);
                assignServer(ch1);
                System.out.println("total load after adding user = "+totalLoad(load));
                mainAlgo();
            }
        }while(ch==1);
    }
    
    //funtion to assign server
    public static void assignServer(int ch){
        
        Random rn = new Random();
        int range = m ;
        while(true){
            int randomNum =  rn.nextInt(range) + 0;
            int s = randomNum;
            System.out.println("s = "+s);
            if(server[s].n < thresh){
                P[ch][s] = 1;
                server[s].n += 1;
                server[s].V.add(ch);
                server[s].L = readTransferLoad(server[s].n, server[s].V, s);
                double r_load = readLoad(server[s].n, server[s].V, s, server[s].L);
                double w_load = writeLoad(server[s].n, s, server[s].L, server[s].nr);
                load[s] = r_load + w_load; 
                break;
            }
        }
    }
    //function to decide relationship of new user
    public static void createRelation(int user)
    {
        rel[user] = new ArrayList<>();
        for(int i = 0; i < n-1; i++)
        {
            
            int random;
            random = (Math.random() <= 0.5)? 0 : 1;
            if(random == 1)
            {
                if(rel[user].size() < (user/10)){
                    rel[user].add(i);
                    rel[i].add(user);
                }else
                    break;
                
            }
        }
    }
    

    private static void printServer(int i) throws IOException {
        bw.write("<--------- Server "+i+" --------->\n");
        System.out.println("<--------- Server "+i+" --------->");
        bw.write("No. of users: "+server[i].n+"\n");
        System.out.println("No. of users: "+server[i].n);
        bw.write("Users List: ");
        System.out.print("Users List: ");
        for(int j=0;j<server[i].n;j++){
            bw.write(server[i].V.get(j)+" ");
            System.out.print(server[i].V.get(j)+" ");
        }
        System.out.println();
        bw.write("\nNo. of replicated users: "+server[i].nr+"\n");
        System.out.println("No. of replicated users: "+server[i].nr);
        bw.write("Replicated Users List: ");
        System.out.print("Replicated Users List: ");
        for(int j=0;j<server[i].nr;j++){
                bw.write(server[i].Vr.get(j)+" ");
                System.out.print(server[i].Vr.get(j)+" ");
        }
        System.out.println();
        bw.write("\nLoad at server: "+load[i]+"\n\n");
        System.out.println("Load at server: "+load[i]);
        System.out.println();
    }

    private static void printResult() throws IOException {
        bw.write("************ RESULT ************\n");
        System.out.println("************ RESULT ************");
        DecimalFormat df = new DecimalFormat("#.####");
        df.setRoundingMode(RoundingMode.CEILING);
        bw.write("Server\t\t");
        System.out.print("Server\t\t");
        for(int i=0;i<m;i++){
            bw.write(i+"\t\t\t");
            System.out.print(i+"\t\t\t");
        }
        bw.write("\nLoad\t\t");
        System.out.print("\nLoad\t\t");
        for(int i=0;i<m;i++){
            bw.write(df.format(load[i])+"\t\t\t");
            System.out.print(df.format(load[i])+"\t\t\t");
        }
        bw.write("\n#users\t\t");
        System.out.print("\n#users\t\t");
        for(int i=0;i<m;i++){
            bw.write(server[i].n+"\t\t\t");
            System.out.print(server[i].n+"\t\t\t");
        }
        bw.write("\nUsers\t\t");
        System.out.print("\nUsers\t\t");
        for(int i=0;i<m;i++){
            for(int j=0;j<server[i].n;j++){
                bw.write(server[i].V.get(j)+" ");
                System.out.print(server[i].V.get(j)+" ");
            }
            bw.write("\t\t\t");
            System.out.print("\t\t\t");
        }
        bw.write("\n#replicated\t");
        System.out.print("\n#replicated\t");
        for(int i=0;i<m;i++){
            bw.write(server[i].nr+"\t\t\t");
            System.out.print(server[i].nr+"\t\t\t");
        }
        bw.write("\nReplicated\t");
        System.out.print("\nReplicated\t");
        for(int i=0;i<m;i++){
            for(int j=0;j<server[i].nr;j++){
                bw.write(server[i].Vr.get(j)+" ");
                System.out.print(server[i].Vr.get(j)+" ");
            }
            bw.write("\t\t");
            System.out.print("\t\t");
        }
        bw.write("\n\n");
        System.out.println();
        double tload = totalLoad(load);
        bw.write("Total load of the system: "+tload);
        System.out.println("Total load of the system: "+tload);
    }
}

class Server{
    public int n;   //No. of users assigned to server
    public List<Integer> V;    //Set of users assigned to server
    public double L[][];   //Load of reading and transferring user's data to other servers
    public int nr;   //No. of users replicated on server
    public List<Integer> Vr;    //Set of users replicated on server
}