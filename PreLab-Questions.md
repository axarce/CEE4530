## Laboratory 1 PreLab Questions
### Anthony Arce

1. Why are contact lenses hazardous in the laboratory?
    Contact lenses are hazardous in the laboratory as they can trap foreign chemicals inside the cornea and prevents effective eye washing.
2. What is the minimum information needed on the label for each chemical? When are Right-To-Know labels required?
  The minimum information provided on the label should be its full chemical name, concentration, and date prepared. Right-To-Know labels are required when a chemical is a hazardous material in a solution, solid, or liquid form in a concentration of greater than one percent or one-tenth of one percent if the material in question is a carcinogen.
3. Why is it important to label a bottle even if it only contains distilled water?
  When it comes to disposal, it makes it much simpler to have it labeled as you or anyone else would be able to follow instructions for disposal regardless if it is distilled water.
4. Find an SDS for sodium nitrate.
    1. Who created the SDS?
      Sciencelab.com, Inc. a chemical and laboratory equipment company.
    2. What is the solubility of sodium nitrate in water?
      Dependent on the temperature. 92.1 g Sodium Nitrate / 100 ml water at 25 degrees C. Increase in solubility at higher temperatures.
    3. Is sodium nitrate carcinogenic?
      No, sodium nitrate is not carcinogenic.
    4. What is the LD50 oral rat?
      The LD50 orally for rats is 1267 mg/kg.
    5. How much sodium nitrate would you have to ingest to give a 50% chance of death (estimate from available LD50 data).
    ```python
    from aide_design.play import *
    myWeight = 160*u.lb
    lb2kg = 0.41*u.kg/u.lb
    ld50 = 1267*u.mg/u.kg
    ld50self = myWeight*lb2kg*ld50
    print(ld50self)
    #insert calculations here
    ```
    6. How much of a 1 M solution would you have to ingest to give a 50% chance of death?
    ```python
    from aide_design.play import *
    ld50self = 83.12*u.g
    molarity = 1*u.mole/u.liter
    NaNO3MolWt = 84.99*u.g/u.mole
    ld50vol = ld50self/(molarity*NaNO3MolWt)
    print(ld50vol)
    #insert calculations here
    ```
    7. Are there any chronic effects of exposure to sodium nitrate?
      Chronic exposure might be toxic to blood and produce organ damage.
5. You are in the laboratory preparing chemical solutions for an experiment and it is lunchtime. You decide to go to CTB to eat. What must you do before leaving the laboratory?
  You must wash your hands prior to leaving the laboratory.
