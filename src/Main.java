import java.io.*;
import java.util.*;

public class Main {
    public static void main(String[] args) throws IOException {
        // params
        List<Object> params = new ArrayList<>();

        // todo add modulo function

        try {
            File myObj = new File("/Users/bartek/intellij-workspace/psm_lab1/src/params.txt");
            Scanner myReader = new Scanner(myObj);
            int lineCount = 0;
            while (myReader.hasNextLine()) {
                String data = myReader.nextLine();
                if (lineCount > 8) {
                    params.add(Double.parseDouble(data));
                } else {
                    params.add(Integer.parseInt(data));
                }
                lineCount++;
            }
            myReader.close();
        } catch (FileNotFoundException exception) {
            System.out.println("An error occurred.");
            exception.printStackTrace();
        }

        final int n = (int) params.get(0);
        final int m = (int) params.get(1);
        final int e = (int) params.get(2);
        final int f = (int) params.get(3);
        final int T_0 = (int) params.get(4);
        final int S_0 = (int) params.get(5);
        final int S_d = (int) params.get(6);
        final int S_out = (int) params.get(7);
        final int S_xyz = (int) params.get(8);
        final double L = (double) params.get(9);
        final double a = (double) params.get(10);
        final double R = (double) params.get(11);
        final double tau = (double) params.get(12);

        final double k = 0.00831;

        // 2.1

        int N = (int) Math.pow(n, 3);

        double[] b0 = new double[] {a, 0, 0};
        double[] b1 = new double[] {a/2, a*Math.sqrt(3)/2, 0};
        double[] b2 = new double[] {a/2, a*Math.sqrt(3)/6, a*Math.sqrt((double) 2/3)};

        List<double[]> coordinates = new ArrayList<>();
        double constValue = (double) (n - 1)/2;
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                for (int l = 0; l < n; l++) {
                    double x = (i - constValue)*b0[0] + (i - constValue)*b1[0] + (i - constValue)*b2[0];
                    double y = (j - constValue)*b0[1] + (j - constValue)*b1[1] + (j - constValue)*b2[1];
                    double z = (l - constValue)*b0[2] + (l - constValue)*b1[2] + (l - constValue)*b2[2];
                    coordinates.add(new double[]{x, y, z});
                }
            }
        }

        System.out.println(Arrays.toString(coordinates.get(2)));

        BufferedWriter writer = new BufferedWriter(new FileWriter("coordinates.xyz"));
        writer.write(1000 + "\n");
        writer.write("comment\n");
        for (double[] coordinate : coordinates) {
            writer.write("Ar" + "\t" + coordinate[0] + "\t" + coordinate[1] + "\t" + coordinate[2] + "\n");
        }
        writer.close();

        List<double[]> energies = new ArrayList<>(N);
        Random random = new Random();
        for (int i = 0; i < N; i++) {
            double energyX = (-1.0/2)*k*T_0*Math.log(random.nextDouble(0, 1));
            double energyY = (-1.0/2)*k*T_0*Math.log(random.nextDouble(0, 1));
            double energyZ = (-1.0/2)*k*T_0*Math.log(random.nextDouble(0, 1));
            energies.add(new double[]{energyX, energyY, energyZ});
        }

        List<double[]> momentums = new ArrayList<>(N);
        for (double[] energy : energies) {
            double pendX = (random.nextDouble(0, 1) < 0.5 ? -1 : 1)*Math.sqrt(2*m*energy[0]);
            double pendY = (random.nextDouble(0, 1) < 0.5 ? -1 : 1)*Math.sqrt(2*m*energy[1]);
            double pendZ = (random.nextDouble(0, 1) < 0.5 ? -1 : 1)*Math.sqrt(2*m*energy[2]);
            momentums.add(new double[]{pendX, pendY, pendZ});
        }

        writer = new BufferedWriter(new FileWriter("energies"));
        for (int i = 0; i < 3; i++) {
            for (double[] momentum : momentums) {
                writer.write(momentum[i] + "\t");
            }
            writer.write("\n");
        }
        writer.close();

        // 2.2

        List<Double> V_Pij = new ArrayList<>(N);
        List<double[]> radius_dif = new ArrayList<>(N);
        double V_P = 0;
        for (int i = 0; i < N; i++) {
            double V_PSum = 0;
            radius_dif.add(new double[N]);
            for (int j = 0; j < N; j++) {
                if (i != j) {
                    double r1 = coordinates.get(i)[0] - coordinates.get(j)[0];
                    double r2 = coordinates.get(i)[1] - coordinates.get(j)[1];
                    double r3 = coordinates.get(i)[2] - coordinates.get(j)[2];
                    double r = Math.sqrt(Math.pow(r1, 2) + Math.pow(r2, 2) + Math.pow(r3, 2));
//                    System.out.println(r);
                    V_PSum += e*(Math.pow(R/r, 12) - 2*Math.pow(R/r, 6));
                    radius_dif.get(i)[j] = r;
                } else {
                    radius_dif.get(i)[j] = 0;
                }
            }
            V_Pij.add(V_PSum);
            V_P += V_PSum;
        }

        System.out.println(V_P);

        List<Double> V_Sij = new ArrayList<>(N);
        List<Double> radius = new ArrayList<>(N);
        double V_S = 0;
        for (int i = 0; i < N; i++) {
            double r = Math.sqrt(Math.pow(coordinates.get(i)[0], 2)
                    + Math.pow(coordinates.get(i)[1], 2)
                    + Math.pow(coordinates.get(i)[2], 2));
            radius.add(r);
            if (r < L) {
                V_Sij.add(0.0);
            } else {
                double V_Sij_val = ((double) 1/2)*f*Math.pow(r - L, 2);
                V_Sij.add(V_Sij_val);
                V_S += V_Sij_val;
            }
        }

        System.out.println(V_S);

        double V = V_P + V_S;

        List<List<double[]>> F_Pij = new ArrayList<>(N);
        double F_P = 0;
        for (int i = 0; i < N; i++) {
            F_Pij.add(new ArrayList<>(N - 1));
            for (int j = 0; j < N; j++) {
                if (i != j) {
                    double r = radius_dif.get(i)[j];
                    double r_x = coordinates.get(i)[0] - coordinates.get(j)[0];
                    double r_y = coordinates.get(i)[1] - coordinates.get(j)[1];
                    double r_z = coordinates.get(i)[2] - coordinates.get(j)[2];

                    double stupidParam = Math.pow(R/r, 12) - 2*Math.pow(R/r, 6);

                    double F_x = 12.0*e*stupidParam*r_x/Math.pow(r, 2);
                    double F_y = 12.0*e*stupidParam*r_y/Math.pow(r, 2);
                    double F_z = 12.0*e*stupidParam*r_z/Math.pow(r, 2);

                    F_Pij.get(i).add(new double[]{F_x, F_y, F_z});
                }
            }
        }

        System.out.println(F_P);

        List<double[]> F_Si = new ArrayList<>(N);
        for (int i = 0; i < N; i++) {
            double r = radius.get(i);
            if (r < L) {
                F_Si.add(new double[]{0, 0, 0});
            } else {
                double F_x = f*(L - r)*coordinates.get(i)[0]/r;
                double F_y = f*(L - r)*coordinates.get(i)[1]/r;
                double F_z = f*(L - r)*coordinates.get(i)[2]/r;
                F_Si.add(new double[]{F_x, F_y, F_z});
            }
        }

        double P = 0;
        for (int i = 0; i < N; i++) {
            P += Math.sqrt(Math.pow(F_Si.get(i)[0], 2) + Math.pow(F_Si.get(i)[1], 2) + Math.pow(F_Si.get(i)[2], 2));
        }
        P  = P/(4*Math.PI*Math.pow(L, 2));
        System.out.println(P);

        // 2.3

        double H = 0;
        for (int i = 0; i < N; i++) {
            H += Math.abs(Math.pow(momentums.get(i)[0], 2)
                    + Math.pow(momentums.get(i)[1], 2)
                    + Math.pow(momentums.get(i)[2], 2));
        }
        H = H/(2*m);
        H += V;

        System.out.println(H);



        for (int s = 0; s < S_0 + S_d; s++) {

        }
        

    }
}
