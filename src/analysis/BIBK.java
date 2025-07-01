package analysis;

import core.Monomer;
import core.Polymer;
import core.Vector;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

public class BIBK {
    public Map<Integer, Double> calcBIBK(List<Polymer> polymers) {
        TreeMap<Integer, Double> dotProductSums = new TreeMap<>();
        TreeMap<Integer, Integer> counts = new TreeMap<>();

        for (Polymer chainState : polymers) {
            ArrayList<Monomer> monomers = chainState.getMonomers();
            if (monomers.size() < 2) {
                continue;
            }

            ArrayList<Vector> tangents = new ArrayList<>();
            for (int i = 0; i < monomers.size() - 1; i++) {
                Vector r_i_iplus1 = monomers.get(i + 1).getPosition().sub(monomers.get(i).getPosition());
                tangents.add(r_i_iplus1.normalize());
            }

            for (int i = 0; i < tangents.size(); i++) {
                for (int j = i; j < tangents.size(); j++) {
                    int separation = j - i;
                    if (separation == 0) {
                        continue;
                    }

                    double dotProduct = tangents.get(i).dot(tangents.get(j));

                    dotProductSums.put(separation, dotProductSums.getOrDefault(separation, 0.0) + dotProduct);
                    counts.put(separation, counts.getOrDefault(separation, 0) + 1);
                }
            }
        }

        TreeMap<Integer, Double> correlationFunction = new TreeMap<>();
        for (int separation : dotProductSums.keySet()) {
            correlationFunction.put(separation, dotProductSums.get(separation) / counts.get(separation));
        }

        return correlationFunction;
    }
}
