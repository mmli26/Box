package analysis;

import core.Vector;
import core.Monomer;
import core.Polymer;

import java.util.List;
import java.util.TreeMap;
import java.util.ArrayList;
import java.util.Map;

public class R2S {
    public Map<Integer, Double> calcR2S(List<Polymer> trajectory){
        TreeMap<Integer, Double> distanceSqSums = new TreeMap<>();
        TreeMap<Integer, Integer> count = new TreeMap<>();

        for (Polymer p : trajectory){
            ArrayList<Monomer> monomers = p.getMonomers();
            int chainLength = monomers.size();
            if (chainLength < 2){
                continue;
            }
            for (int s = 1; s < chainLength; s++){
                for (int i = 0; i < chainLength - s; i++){
                    int j = i + s;

                    Monomer m1 = monomers.get(i);
                    Monomer m2 = monomers.get(j);

                    Vector r_ij = m2.getPosition().sub(m1.getPosition());
                    double distSq = r_ij.magnitude() * r_ij.magnitude();

                    distanceSqSums.put(s, distanceSqSums.getOrDefault(s, 0.0) + distSq);
                    count.put(s, count.getOrDefault(s, 0) + 1);
                }
            }
        }
        TreeMap<Integer, Double> r2s = new TreeMap<>();
        for (int s : distanceSqSums.keySet()){
            r2s.put(s, distanceSqSums.get(s) / count.get(s));
        }
        return r2s;
    }
}
