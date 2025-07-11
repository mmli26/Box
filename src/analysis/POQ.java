package analysis;

import core.Monomer;
import core.Polymer;
import core.Vector;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.concurrent.*;

public class POQ {
    public Map<Double, Double> calcPOQ(List<Polymer> trajectory, List<Double> qValues) {
        ConcurrentMap<Double, Double> pqSums = new ConcurrentHashMap<>();
        for (double q : qValues) {
            pqSums.put(q, 0.0);
        }

        int nThreads = Runtime.getRuntime().availableProcessors();
        ExecutorService executor = Executors.newFixedThreadPool(nThreads);
        List<Future<?>> futures = new CopyOnWriteArrayList<>();

        int chunkSize = (int) Math.ceil((double) trajectory.size() / nThreads);

        for (int i = 0; i < nThreads; i++) {
            final int start = i * chunkSize;
            final int end = Math.min(start + chunkSize, trajectory.size());

            if (start >= end) continue;

            List<Polymer> subTrajectory = trajectory.subList(start, end);

            Future<?> future = executor.submit(() -> {
                for (Polymer snapshot : subTrajectory) {
                    ArrayList<Monomer> monomers = snapshot.getMonomers();
                    int chainLength = monomers.size();
                    if (chainLength == 0) continue;

                    for (int m1 = 0; m1 < chainLength; m1++) {
                        Monomer monomer_i = monomers.get(m1);
                        for (int m2 = 0; m2 < chainLength; m2++) {
                            Monomer monomer_j = monomers.get(m2);

                            double r_ij;
                            if (m1 == m2) {
                                r_ij = 0.0;
                            } else {
                                Vector displacement = monomer_i.getPosition().sub(monomer_j.getPosition());
                                r_ij = displacement.magnitude();
                            }

                            for (double q : qValues) {
                                double term;
                                if (r_ij < 1e-9) {
                                    term = 1.0;
                                } else {
                                    term = Math.sin(q * r_ij) / (q * r_ij);
                                }
                                pqSums.compute(q, (key, val) -> (val == null) ? term : val + term);
                            }
                        }
                    }
                }
            });
            futures.add(future);
        }

        for (Future<?> f : futures) {
            try {
                f.get();
            } catch (InterruptedException | ExecutionException e) {
                e.printStackTrace();
            }
        }
        executor.shutdown();

        Map<Double, Double> finalPq = new TreeMap<>();
        int N = trajectory.get(0).getMonomers().size();
        double normalizationFactor = (double) N * N * trajectory.size();
        if (normalizationFactor == 0) return finalPq;

        for (double q : pqSums.keySet()) {
            finalPq.put(q, pqSums.get(q) / normalizationFactor);
        }

        return finalPq;
    }
}