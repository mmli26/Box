package analysis;

import core.Polymer;
import core.Vector;

import java.util.List;
import java.util.Map;
import java.util.TreeMap;

public class MSD_MON {
    public Map<Integer, Double> calcMsd(List<Polymer> trajectory) {
        TreeMap<Integer, Double> msdSums = new TreeMap<>();
        TreeMap<Integer, Integer> counts = new TreeMap<>();

        int nSnapshots = trajectory.size();
        if (nSnapshots < 2) {
            return new TreeMap<>();
        }

        int nParticles = trajectory.get(0).getMonomers().size();
        if (nParticles == 0) {
            return new TreeMap<>();
        }

        for (int dt = 1; dt < nSnapshots; dt++) {
            double currentMsdSum = 0.0;
            int currentCount = 0;

            for (int i = 0; i < nSnapshots - dt; i++) {
                int j = i + dt;

                Polymer snapshot_i = trajectory.get(i);
                Polymer snapshot_j = trajectory.get(j);

                for (int p = 0; p < nParticles; p++) {
                    Vector pos_i = snapshot_i.getMonomers().get(p).getPosition();
                    Vector pos_j = snapshot_j.getMonomers().get(p).getPosition();

                    Vector displacement = pos_j.sub(pos_i);
                    currentMsdSum += displacement.dot(displacement);
                    currentCount++;
                }
            }

            if (currentCount > 0) {
                msdSums.put(dt, currentMsdSum);
                counts.put(dt, currentCount);
            }
        }

        TreeMap<Integer, Double> finalMsd = new TreeMap<>();
        for (int dt : msdSums.keySet()) {
            finalMsd.put(dt, msdSums.get(dt) / counts.get(dt));
        }

        return finalMsd;
    }
}
