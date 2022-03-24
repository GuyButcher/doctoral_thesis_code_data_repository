% Generates the properties variable that contains the characteristics of
% each amino acid and their sidechains.
%
% KD Hydropathy source article:
% "https://www.sciencedirect.com/science/article/abs/pii/0022283682905150"
% Value is positive for hydrophobic regions and negative for hydrophilic
% regions.
% Additional information taken from : "https://home.hiroshima-u.ac.jp/kei/IdentityX/picts/BE-hydrophobicity.pdf"
%
% Charges source book: "https://www.ncbi.nlm.nih.gov/books/NBK9879/"
%
% actually taken from 
% "https://en.wikipedia.org/wiki/Amino_acid#Physicochemical_properties_of_amino_acids"
% and "https://en.wikipedia.org/wiki/Proteinogenic_amino_acid"
% and "https://en.wikipedia.org/wiki/Hydrophilicity_plot"

residue_properties.name_long = strings(20,1);
residue_properties.name_short = strings(20,1);
residue_properties.name_letter = strings(20,1);
residue_properties.mass = zeros(20,1);
residue_properties.charge = zeros(20,1);
residue_properties.hydrophobicity_kd = zeros(20,1);

residue_properties.name_long(1)=            "Alanine"; 
residue_properties.name_short(1) =          "Ala";
residue_properties.name_letter(1) =         "A";
residue_properties.mass(1) =                0;
residue_properties.charge(1) =              0;
residue_properties.hydrophobicity_kd(1) =   1.8;

residue_properties.name_long(2)=            "Arginine";
residue_properties.name_short(2) =          "Arg";
residue_properties.name_letter(2) =         "R";
residue_properties.mass(2) =                0;
residue_properties.charge(2) =              1;
residue_properties.hydrophobicity_kd(2) =   -4.5;

residue_properties.name_long(3)=            "Asparagine";
residue_properties.name_short(3) =          "Asn";
residue_properties.name_letter(3) =         "N";
residue_properties.mass(3) =                0;
residue_properties.charge(3) =              0;
residue_properties.hydrophobicity_kd(3) =   -3.5;

residue_properties.name_long(4)=            "Aspartic Acid";
residue_properties.name_short(4) =          "Asp";
residue_properties.name_letter(4) =         "D";
residue_properties.mass(4) =                0;
residue_properties.charge(4) =              -1;
residue_properties.hydrophobicity_kd(4) =   -3.5;

residue_properties.name_long(5)=            "Cysteine";
residue_properties.name_short(5) =          "Cys";
residue_properties.name_letter(5) =         "C";
residue_properties.mass(5) =                0;
residue_properties.charge(5) =              0;
residue_properties.hydrophobicity_kd(5) =   2.5;

residue_properties.name_long(6)=            "Glutamine";
residue_properties.name_short(6) =          "Gln";
residue_properties.name_letter(6) =         "Q";
residue_properties.mass(6) =                0;
residue_properties.charge(6) =              0;
residue_properties.hydrophobicity_kd(6) =   -3.5;

residue_properties.name_long(7)=            "Glutamic Acid";
residue_properties.name_short(7) =          "Glu";
residue_properties.name_letter(7) =         "E";
residue_properties.mass(7) =                0;
residue_properties.charge(7) =              -1;
residue_properties.hydrophobicity_kd(7) =   -3.5;

residue_properties.name_long(8)=            "Glycine";
residue_properties.name_short(8) =          "Gly";
residue_properties.name_letter(8) =         "G";
residue_properties.mass(8) =                0;
residue_properties.charge(8) =              0;
residue_properties.hydrophobicity_kd(8) =   -0.4;

residue_properties.name_long(9)=            "Histidine";
residue_properties.name_short(9) =          "His";
residue_properties.name_letter(9) =         "H";
residue_properties.mass(9) =                0;
residue_properties.charge(9) =              1;
residue_properties.hydrophobicity_kd(9) =   -3.2;

residue_properties.name_long(10)=            "Isoleucine";
residue_properties.name_short(10) =          "Ile";
residue_properties.name_letter(10) =         "I";
residue_properties.mass(10) =                0;
residue_properties.charge(10) =              0;
residue_properties.hydrophobicity_kd(10) =   4.5;

residue_properties.name_long(11)=            "Leucine";
residue_properties.name_short(11) =          "Leu";
residue_properties.name_letter(11) =         "L";
residue_properties.mass(11) =                0;
residue_properties.charge(11) =              0;
residue_properties.hydrophobicity_kd(11) =   3.8;

residue_properties.name_long(12)=            "Lysine";
residue_properties.name_short(12) =          "Lys";
residue_properties.name_letter(12) =         "K";
residue_properties.mass(12) =                0;
residue_properties.charge(12) =              1;
residue_properties.hydrophobicity_kd(12) =   -3.9;

residue_properties.name_long(13)=            "Methionine";
residue_properties.name_short(13) =          "Met";
residue_properties.name_letter(13) =         "M";
residue_properties.mass(13) =                0;
residue_properties.charge(13) =              0;
residue_properties.hydrophobicity_kd(13) =   1.9;

residue_properties.name_long(14)=            "Phenylalanine";
residue_properties.name_short(14) =          "Phe";
residue_properties.name_letter(14) =         "F";
residue_properties.mass(14) =                0;
residue_properties.charge(14) =              0;
residue_properties.hydrophobicity_kd(14) =   2.8;

residue_properties.name_long(15)=            "Proline";
residue_properties.name_short(15) =          "Pro";
residue_properties.name_letter(15) =         "P";
residue_properties.mass(15) =                0;
residue_properties.charge(15) =              0;
residue_properties.hydrophobicity_kd(15) =   -1.6;

residue_properties.name_long(16)=            "Serine";
residue_properties.name_short(16) =          "Ser";
residue_properties.name_letter(16) =         "S";
residue_properties.mass(16) =                0;
residue_properties.charge(16) =              0;
residue_properties.hydrophobicity_kd(16) =   -0.8;

residue_properties.name_long(17)=            "Threonine";
residue_properties.name_short(17) =          "Thr";
residue_properties.name_letter(17) =         "T";
residue_properties.mass(17) =                0;
residue_properties.charge(17) =              0;
residue_properties.hydrophobicity_kd(17) =   -0.7;

residue_properties.name_long(18)=            "Tryptophan";
residue_properties.name_short(18) =          "Trp";
residue_properties.name_letter(18) =         "W";
residue_properties.mass(18) =                0;
residue_properties.charge(18) =              0;
residue_properties.hydrophobicity_kd(18) =   -0.9;

residue_properties.name_long(19)=            "Tyrosine";
residue_properties.name_short(19) =          "Tyr";
residue_properties.name_letter(19) =         "Y";
residue_properties.mass(19) =                0;
residue_properties.charge(19) =              0;
residue_properties.hydrophobicity_kd(19) =   -1.3;

residue_properties.name_long(20)=            "Valine";
residue_properties.name_short(20) =          "Val";
residue_properties.name_letter(20) =         "V";
residue_properties.mass(20) =                0;
residue_properties.charge(20) =              0;
residue_properties.hydrophobicity_kd(20) =   4.2;