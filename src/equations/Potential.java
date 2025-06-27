package equations;

import core.Monomer;

import java.util.ArrayList;


public interface Potential {
    double calcEnergy(double value);

    String printEnergy(double value);

    void applyForces(ArrayList<Monomer> monomers);
}
