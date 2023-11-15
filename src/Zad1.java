import java.io.*;
import java.util.*;

public class Zad1 {
    public static void main(String[] args) throws IOException {
        List<Object> params = new ArrayList<>();

        File myObj = new File("params.txt");
        Scanner myReader = new Scanner(myObj);
        int lineCount = 0;
        while (myReader.hasNextLine()) {
            String data = myReader.nextLine().split("\\s+")[0];
            if (lineCount > 8) {
                params.add(Double.parseDouble(data));
            } else {
                params.add(Integer.parseInt(data));
            }
            lineCount++;
        }
        myReader.close();

        final int n = (int) params.get(0);
        final int m = (int) params.get(1);
        final int e = (int) params.get(2);
        final int f = (int) params.get(3);
        final int T_0 = (int) params.get(4);
        final int S_0 = (int) params.get(5);
        final int S_d = (int) params.get(6);
        final int S_out = (int) params.get(7);
        final int S_xyz = (int) params.get(8);
        final double a = (double) params.get(10);
        final double R = (double) params.get(11);
        final double tau = (double) params.get(12);

        final double L = 1.23*a*(n-1);
        final double k = 0.00831;
        double time = 0;

        int N = (int) Math.pow(n, 3);

        double[] b0 = new double[] {a, 0, 0};
        double[] b1 = new double[] {a/2, a*Math.sqrt(3)/2, 0};
        double[] b2 = new double[] {a/2, a*Math.sqrt(3)/6, a*Math.sqrt((double) 2/3)};

        double[][] coordinates = new double[N][3];
        double constValue = (double) (n - 1)/2;
        for (int i_0 = 0; i_0 < n; i_0++) {
            for (int i_1 = 0; i_1 < n; i_1++) {
                for (int i_2 = 0; i_2 < n; i_2++) {
                    int i = i_0 + i_1*n + i_2*n*n;
                    getCoordinates(coordinates[i], b0, b1, b2, constValue, i_0, i_1, i_2);
                }
            }
        }

        BufferedWriter writer = new BufferedWriter(new FileWriter("coordinates.xyz"));
        writer.write(N + "\n\n");
        for (double[] coordinate : coordinates) {
            writer.write("Ar" + "\t" + coordinate[0] + "\t" + coordinate[1] + "\t" + coordinate[2] + "\n");
        }
        writer.close();

        double[][] energies = new double[N][3];
        Random random = new Random();
        constValue = -0.5*k*T_0;
        for (double[] energy : energies) {
            energy[0] = constValue*Math.log(random.nextDouble(0, 1));
            energy[1] = constValue*Math.log(random.nextDouble(0, 1));
            energy[2] = constValue*Math.log(random.nextDouble(0, 1));
        }

        double[][] momenta = new double[N][3];
        for (int i = 0; i < N; i++) {
            momenta[i][0] = (random.nextDouble(0, 1) < 0.5 ? -1 : 1)*Math.sqrt(2*m*energies[i][0]);
            momenta[i][1] = (random.nextDouble(0, 1) < 0.5 ? -1 : 1)*Math.sqrt(2*m*energies[i][1]);
            momenta[i][2] = (random.nextDouble(0, 1) < 0.5 ? -1 : 1)*Math.sqrt(2*m*energies[i][2]);
        }

        double momentaShift = sumMoments(momenta)/N;
        shiftMomenta(momenta, momentaShift);

        writer = new BufferedWriter(new FileWriter("energies"));
        for (int i = 0; i < 3; i++) {
            for (double[] momentum : momenta) {
                writer.write(momentum[i] + "\t");
            }
            writer.write("\n");
        }
        writer.close();

        long timeStart = System.currentTimeMillis();
        double[][] F = new double[N][3];
        double[] VAndP = countPotentialsAndForces(N, coordinates, L, f, R, e, F);
        double V = VAndP[0];

        double avrPressure = 0;
        double avrTemperature = 0;
        double avrHamiltonian = 0;

        BufferedWriter argonWriter = new BufferedWriter(new FileWriter("argon.txt"));
        BufferedWriter parametersWriter = new BufferedWriter(new FileWriter("parameters.xyz"));
        BufferedWriter energyWriter = new BufferedWriter(new FileWriter("absolute_energy.xyz"));

        for (int i = 0; i < S_d; i++) {
            double[] potentialAndPressure = equationsOfMotion(coordinates, momenta, tau, F, N, m, L, f, R, e);
            double[] kinEnergies = updateEnergies(N, m, momenta);
            double T = countTemperature(N, k, kinEnergies);
            double H = countHamiltonian(m, V, momenta);
            time += tau;

            if (i > S_0) {
                avrPressure += potentialAndPressure[1];
                avrTemperature += T;
                avrHamiltonian += H;
            }

            if (i%S_out == 0) {
                saveParameters(parametersWriter, time, potentialAndPressure[0], potentialAndPressure[1], T, H);
            }

            if (i%S_xyz == 0) {
                saveCoordinatesAndEnergies(argonWriter, N, coordinates, kinEnergies);
            }

            if (i%10 == 0) {
                saveEnergy(energyWriter, kinEnergies, potentialAndPressure[0]);
            }

            if (i%1000 == 0 && i != 0) {
                System.out.println("\nExecution progress: " + i/100 + "%");
                System.out.println("Absolute energy: " + (Arrays.stream(kinEnergies).sum() + potentialAndPressure[0]));
                long timeStop = System.currentTimeMillis();
                System.out.println("Time: " + (timeStop - timeStart)/1000);
            }
        }

        argonWriter.close();
        parametersWriter.close();
        energyWriter.close();

        long timeStop = System.currentTimeMillis();
        System.out.println("\nSimulation finished in: " + (timeStop - timeStart)/1000);

        System.out.println("\nAverage pressure: " + avrPressure/S_d);
        System.out.println("Average temperature: " + avrTemperature/S_d);
        System.out.println("Average hamiltonian: " + avrHamiltonian/S_d);
    }

    private static void getCoordinates(double[] coordinate, double[] b0, double[] b1, double[] b2, double constValue,
                                       int i0, int i1, int i2) {
        double param0 = i0 - constValue;
        double param1 = i1 - constValue;
        double param2 = i2 - constValue;

        coordinate[0] = param0*b0[0] + param1*b1[0] + param2*b2[0];
        coordinate[1] = param0*b0[1] + param1*b1[1] + param2*b2[1];
        coordinate[2] = param0*b0[2] + param1*b1[2] + param2*b2[2];
    }

    private static double sumMoments(double[][] momenta) {
        return Arrays.stream(momenta)
                .map(momentum -> Arrays.stream(momentum).sum())
                .mapToDouble(Double::doubleValue)
                .sum();
    }

    private static void shiftMomenta(double[][] momenta, double momentaShift) {
        for (double[] momentum : momenta) {
            momentum[0] = momentum[0] - momentaShift;
            momentum[1] = momentum[1] - momentaShift;
            momentum[2] = momentum[2] - momentaShift;
        }
    }

    private static double[] countPotentialsAndForces(int N, double[][] coordinates, double L, int f, double R, int e,
                                                 double[][] F) {
        double[] Vs = new double[N];
        double[] Vp = new double[N];
        double[][] Fp = new double[N][3];
        double[][] Fs = new double[N][3];

        for (int i = 0; i < N; i++) {
            double abs = abs(coordinates[i]);

            if (abs > L) {
                Vs[i] = (f*(abs - L)*(abs - L))/2;
                fillFs(Fs[i], L, f, abs, coordinates[i]);
            }

            for (int j = 0; j < i; j++) {
                abs = abs(subtractVectors(coordinates[i], coordinates[j]));
                double sup = R/abs;
                double constant12 = sup*sup*sup*sup*sup*sup*sup*sup*sup*sup*sup*sup;
                double constant6 = sup*sup*sup*sup*sup*sup;

                Vp[i] += e*(constant12 - 2*constant6);
                fillFp(Fp[i], e, abs, constant12, constant6, coordinates[i], coordinates[j], true);
                fillFp(Fp[j], e, abs, constant12, constant6, coordinates[i], coordinates[j], false);
            }
        }

        sumForces(F, Fs, Fp);
        double pressure = countPressure(Fs, L);
        return new double[]{Arrays.stream(Vs).sum() + Arrays.stream(Vp).sum(), pressure};
    }

    private static void fillFs(double[] Fs, double L, int f, double abs, double[] vector) {
        double constant = f*(L - abs)/abs;
        for (int i = 0; i < 3; i++) {
            Fs[i] = vector[i]*constant;
        }
    }

    private static double abs(double[] vector) {
        return Math.sqrt(vector[0]*vector[0] + vector[1]*vector[1] + vector[2]*vector[2]);
    }

    private static double[] subtractVectors(double[] vector1, double[] vector2) {
        return new double[]{vector1[0] - vector2[0], vector1[1] - vector2[1], vector1[2] - vector2[2]};
    }

    private static void fillFp(double[] Fp, int e, double abs, double constant12, double constant6,
                               double[] vector1, double[] vector2, boolean add) {

        double constant = 12*e*(constant12 - constant6)/(abs*abs);
        double[] vector = subtractVectors(vector1, vector2);

        for (int i = 0; i < 3; i++) {
            double tmp = vector[i]*constant;
            if (add) {
                Fp[i] += tmp;
            } else {
                Fp[i] -= tmp;
            }
        }
    }

    private static void sumForces(double[][] F, double[][] Fs, double[][] Fp) {
        for (int i = 0; i < Fs.length; i++) {
            for (int j = 0; j < 3; j++) {
                F[i][j] = Fs[i][j] + Fp[i][j];
            }
        }
    }

    private static double countPressure(double[][] Fs, double L) {
        return Arrays.stream(Fs)
                .map(Zad1::abs)
                .mapToDouble(Double::doubleValue)
                .sum()/(4*Math.PI*L*L);
    }

    private static double[] equationsOfMotion(double[][] coordinates, double[][] momenta, double tau,
                                          double[][] F, int N, int m, double L, int f, double R, int e) {
        updateMomenta(momenta, F, tau, N);
        updateCoordinate(coordinates, momenta, tau, N, m);
        double[] potentialAndPressure = countPotentialsAndForces(N, coordinates, L, f, R, e, F);
        updateMomenta(momenta, F, tau, N);
        return new double[]{potentialAndPressure[0], potentialAndPressure[1]};
    }

    private static void updateMomenta(double[][] momenta, double[][] F, double tau, int N) {
        double constant = tau/2;

        for (int i = 0; i < N; i++) {
            momenta[i][0] += constant*F[i][0];
            momenta[i][1] += constant*F[i][1];
            momenta[i][2] += constant*F[i][2];
        }
    }

    private static void updateCoordinate(double[][] coordinates, double[][] momenta, double tau, int N, int m) {
        double constant = tau/m;

        for (int i = 0; i < N; i++) {
            coordinates[i][0] += constant*momenta[i][0];
            coordinates[i][1] += constant*momenta[i][1];
            coordinates[i][2] += constant*momenta[i][2];
        }
    }

    private static double[] updateEnergies(int N, double m, double[][] momenta) {
        double[] energies = new double[N];
        double constant = 2*m;

        for (int i = 0; i < N; i++) {
            double norm = abs(momenta[i]);
            energies[i] = norm*norm/constant;
        }

        return energies;
    }

    private static double countTemperature(int N, double k, double[] kinEnergies) {
        return 2*Arrays.stream(kinEnergies).sum()/(3*N*k);
    }

    private static double countHamiltonian(int m, double V, double[][] momenta) {
        return Arrays.stream(momenta)
                .map(Zad1::abs)
                .mapToDouble(Double::doubleValue)
                .map(val -> val*val)
                .sum()/(2*m) + V;
    }

    private static void saveCoordinatesAndEnergies(BufferedWriter writer, int N, double[][] coordinates,
                                                   double[] energies) throws IOException {
        writer.write(N + "\n\n");
        for (int i = 0; i < N; i++) {
            writer.write("Ar" + "\t" + coordinates[i][0] + "\t" + coordinates[i][1] + "\t" + coordinates[i][2] +
                    "\t" + energies[i] + "\n");
        }
    }

    private static void saveParameters(BufferedWriter writer, double time, double V, double P, double T, double H)
            throws IOException {
        writer.write(time + "\t" + V + "\t" + P + "\t" + T + "\t" + H + "\n");
    }

    private static void saveEnergy(BufferedWriter writer, double[] energies, double V) throws IOException {
        double kinEnergiesSum = 0;
        for (double energy : energies) {
            kinEnergiesSum += energy;
        }
        double energiesSum = kinEnergiesSum + V;
        writer.write(kinEnergiesSum + "\t" + V + "\t" + energiesSum + "\n");
    }
}
