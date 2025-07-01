package analysis;

import core.Monomer;
import core.Polymer;
import core.Vector;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

public class POQ {
    public Map<Double, Double> calcPOQ(List<Polymer> trajectory, List<Double> qValues) {
        TreeMap<Double, Double> pqSums = new TreeMap<>();
        for (double q : qValues) {
            pqSums.put(q, 0.0);
        }

        double totalChainCount = 0.0;

        for (Polymer snapshot : trajectory) {
            ArrayList<Monomer> monomers = snapshot.getMonomers();
            int chainLength = monomers.size();
            if (chainLength == 0) {
                continue;
            }
            totalChainCount++;

            for (int i = 0; i < chainLength; i++) {
                Monomer monomer_i = monomers.get(i);

                for (int j = 0; j < chainLength; j++) {
                    Monomer monomer_j = monomers.get(j);

                    double r_ij;
                    if (i == j) {
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
                        pqSums.put(q, pqSums.get(q) + term);
                    }
                }
            }
        }

        Map<Double, Double> finalPq = new TreeMap<>();
        if (totalChainCount == 0) {
            return finalPq;
        }

        int N = trajectory.get(0).getMonomers().size();
        double normalizationFactor = (double) N * N * totalChainCount;

        for (double q : pqSums.keySet()) {
            double sum = pqSums.get(q);
            finalPq.put(q, sum / normalizationFactor);
        }

        return finalPq;
    }
}
