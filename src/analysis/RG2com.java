package analysis;

import core.Monomer;
import core.Polymer;
import core.Vector;

import java.util.ArrayList;
import java.util.List;

public class RG2com {
    public double calcRG2com(List<Polymer> trajectory){
        double rg2sum = 0.0;
        int count = 0;

        for (Polymer c : trajectory){
            ArrayList<Polymer> chains = new ArrayList<>();
            chains.add(c);
            for (Polymer p : chains){
                ArrayList<Monomer> monomers = p.getMonomers();
                Vector v = Vector.ZERO;
                for (Monomer m : monomers){
                    v = v.add(m.getPosition());
                }
                int s = monomers.size();
                Vector com = v.scale(1.0 / s);

                double currRg2Sum = 0.0;
                for (Monomer m : monomers){
                    Vector displacement = m.getPosition().sub(com);
                    currRg2Sum += displacement.dot(displacement);
                }
                rg2sum += currRg2Sum;
                count++;
            }
        }
        if (count ==0){
            return 0.0;
        }
        return rg2sum / count;
    }
}
