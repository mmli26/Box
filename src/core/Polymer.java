package core;

import java.util.ArrayList;

public class Polymer {

    private final int id;
    private final ArrayList<Monomer> monomers;
    private final ArrayList<Bond> bonds;
    private final ArrayList<Angle> angles;

    public Polymer(int id){
        this.id = id;
        this.monomers = new ArrayList<>();
        this.bonds = new ArrayList<>();
        this.angles = new ArrayList<>();
    }

    public void addMonomer(Monomer m){
        this.monomers.add(m);
        int currSize = this.monomers.size();
        if (currSize >= 2) {
            Monomer previousMonomer = this.monomers.get(currSize - 2);
            Bond newBond = new Bond(previousMonomer, m);
            this.bonds.add(newBond);
        }
        if (currSize >= 3) {
            Monomer monomer1 = this.monomers.get(currSize - 3);
            Monomer monomer2 = this.monomers.get(currSize - 2);
            Monomer monomer3 = this.monomers.get(currSize - 1);
            Angle newAngle = new Angle(monomer1, monomer2, monomer3);
            this.angles.add(newAngle);
        }
    }

    public int getId(){
        return id;
    }

    public ArrayList<Monomer> getMonomers(){
        return monomers;
    }

    public ArrayList<Bond> getBonds(){
        return bonds;
    }

    public ArrayList<Angle> getAngles(){
        return angles;
    }
}
