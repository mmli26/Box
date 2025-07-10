import analysis.*;
import core.Monomer;
import core.Polymer;
import core.Vector;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Scanner;
import java.util.TreeMap;

public class Main {

    public static void main(String[] args) {
        Scanner scanner = new Scanner(System.in);
        Map<String, String> params = new TreeMap<>();
        System.out.println("Please enter simulation parameters (e.g., 'NN 40', 'directory ./coord/'):");
        System.out.println("Enter 'RUN' when you are finished.");

        while (true) {
            String line = scanner.nextLine();
            if (line.equalsIgnoreCase("RUN")) {
                break;
            }
            String[] parts = line.split("\\s+");
            if (parts.length >= 2) {
                String value = parts[1];
                for (int i = 2; i < parts.length; i++) {
                    value += " " + parts[i];
                }
                params.put(parts[0], value);
            }
        }

        String directory = params.getOrDefault("directory", "./");
        int nn = Integer.parseInt(params.getOrDefault("NN", "0"));
        int initStep = Integer.parseInt(params.getOrDefault("initStep", "0"));
        int endStep = Integer.parseInt(params.getOrDefault("endStep", "0"));
        int dumpStep = Integer.parseInt(params.getOrDefault("dumpStep", "1"));
        String analysisType = params.getOrDefault("analysis", "all");

        System.out.println("Loading trajectory...");
        List<Polymer> trajectory = loadTrajectory(directory, initStep, endStep, dumpStep, nn);
        if (trajectory.isEmpty()) {
            System.err.println("Failed to load any valid trajectory data. Please check your parameters and file format. Exiting.");
            return;
        }

        System.out.println("Successfully loaded " + trajectory.size() + " polymer snapshots.");

        System.out.println("Running analysis: " + analysisType.toUpperCase());
        runAnalysis(analysisType, trajectory, nn);

        System.out.println("Analysis complete. Output written to " + analysisType.toLowerCase() + ".txt");
    }

    /**
     * Reads all specified dump files and constructs a list of polymer states (trajectory).
     */
    public static List<Polymer> loadTrajectory(String directory, int initStep, int endStep, int dumpStep, int nn) {
        List<Polymer> trajectory = new ArrayList<>();
        for (int step = initStep; step <= endStep; step += dumpStep) {
            String filePath = String.format("%s/dump.%09d.txt", directory, step);
            try {
                List<Polymer> snapshot = readDumpFile(filePath, nn);
                if (!snapshot.isEmpty()) {
                    trajectory.addAll(snapshot);
                }
            } catch (IOException e) {
                System.err.println("Warning: Could not find or read file " + filePath + ". Skipping.");
            }
        }
        return trajectory;
    }

    /**
     * Parses a single LAMMPS dump file. This version is robust and specifically
     * handles the 9-column format (id mol type x y z xu yu zu).
     */
    public static List<Polymer> readDumpFile(String filePath, int nn) throws IOException {
        List<Monomer> monomers = new ArrayList<>();
        try (BufferedReader reader = new BufferedReader(new FileReader(filePath))) {
            String line;
            int atomCount = 0;

            while ((line = reader.readLine()) != null && !line.trim().equals("ITEM: NUMBER OF ATOMS")) {}
            if ((line = reader.readLine()) != null) {
                atomCount = Integer.parseInt(line.trim());
            }

            while ((line = reader.readLine()) != null && !line.startsWith("ITEM: ATOMS")) {}

            for (int i = 0; i < atomCount; i++) {
                line = reader.readLine();
                if (line == null || line.trim().isEmpty()) {
                    continue;
                }

                String[] parts = line.trim().split("\\s+");

                if (parts.length < 9) {
                    System.err.println("Warning: Skipping malformed line in " + filePath + ": " + line);
                    continue;
                }

                int id = Integer.parseInt(parts[0]);
                int type = Integer.parseInt(parts[2]);

                double x = Double.parseDouble(parts[6]);
                double y = Double.parseDouble(parts[7]);
                double z = Double.parseDouble(parts[8]);

                monomers.add(new Monomer(new Vector(x, y, z), Vector.ZERO, Vector.ZERO, 1.0, type, id));
            }
        }

        monomers.sort((m1, m2) -> Integer.compare(m1.getIdx(), m2.getIdx()));

        List<Polymer> polymers = new ArrayList<>();
        if (nn > 0 && !monomers.isEmpty()) {
            int nChains = monomers.size() / nn;
            for (int i = 0; i < nChains; i++) {
                Polymer polymer = new Polymer(i + 1);
                for (int j = 0; j < nn; j++) {
                    int monomerIndex = i * nn + j;
                    if (monomerIndex < monomers.size()) {
                        polymer.addMonomer(monomers.get(monomerIndex));
                    }
                }
                polymers.add(polymer);
            }
        }
        return polymers;
    }

    /**
     * Runs the selected analysis and writes the output to a file.
     */
    public static void runAnalysis(String analysisType, List<Polymer> trajectory, int nn) {
        String outputFileName = analysisType.toLowerCase() + ".txt";
        try (FileWriter writer = new FileWriter(outputFileName)) {
            switch (analysisType.toUpperCase()) {
                case "RG2COM":
                    RG2com rg2com = new RG2com();
                    double resultRg2 = rg2com.calcRG2com(trajectory);
                    writer.write("#<Rg^2> " + resultRg2 + "\n");
                    break;
                case "R2S":
                    R2S r2s = new R2S();
                    Map<Integer, Double> resultR2S = r2s.calcR2S(trajectory);
                    for (Map.Entry<Integer, Double> entry : resultR2S.entrySet()) {
                        writer.write(entry.getKey() + " " + entry.getValue() + "\n");
                    }
                    break;
                case "BIBK":
                    BIBK bibk = new BIBK();
                    Map<Integer, Double> resultBibk = bibk.calcBIBK(trajectory);
                    for (Map.Entry<Integer, Double> entry : resultBibk.entrySet()) {
                        writer.write(entry.getKey() + " " + entry.getValue() + "\n");
                    }
                    break;
                case "MSD":
                    MSD_MON msd = new MSD_MON();
                    Map<Integer, Double> resultMsd = msd.calcMsd(trajectory);
                    for (Map.Entry<Integer, Double> entry : resultMsd.entrySet()) {
                        writer.write(entry.getKey() + " " + entry.getValue() + "\n");
                    }
                    break;
                case "POQ":
                    POQ poq = new POQ();
                    List<Double> qValues = new ArrayList<>();
                    for (int i = 0; i < 65; i++) {
                        qValues.add(0.01 * Math.exp((double)i / 10.0));
                    }
                    Map<Double, Double> resultPoq = poq.calcPOQ(trajectory, qValues);
                    for (Map.Entry<Double, Double> entry : resultPoq.entrySet()) {
                        writer.write(entry.getKey() + " " + entry.getValue() + "\n");
                    }
                    break;
                default:
                    System.err.println("Unknown analysis type: " + analysisType);
                    System.err.println("Available types: RG2COM, R2S, BIBK, MSD, POQ");
                    return;
            }
        } catch (IOException e) {
            System.err.println("Error writing to output file " + outputFileName);
            e.printStackTrace();
        }
    }
}