#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define MAX_ATOMS 10000

typedef struct {
  int atom_no;
  char atom_name[5];
  char res_name[4];
  int res_no;
  double x, y, z;
} Atom;

// Sidechain filter logic
int is_sidechain(const char *res, const char *name) {
  if (strcmp(res, "GLY") == 0) {
    return (strcmp(name, "CA") == 0); // Glycine ke liye sirf CA
  }
  if (strcmp(name, "N") == 0 || strcmp(name, "CA") == 0 || 
      strcmp(name, "C") == 0 || strcmp(name, "O") == 0) {
    return 0; // Backbone atoms ko exclude kar rahe hain
  }
  return 1;
}

int main() {
  FILE *infile = fopen("C:/Users/nandini saha/Downloads/Vax D lab/3N42_no_ANISOU with Fchain only.pdb", "r");
  FILE *outfile = fopen("sidechain_matrix.csv", "w");
  
  if (!infile) {
    printf("Error: PDB file nahi mili!\n");
    return 1;
  }
  
  Atom *atoms = malloc(MAX_ATOMS * sizeof(Atom));
  char line[100];
  int count = 0;
  
  // 1. PDB File Read karna aur Sidechain atoms store karna
  while (fgets(line, sizeof(line), infile)) {
    if (strncmp(line, "ATOM", 4) == 0) {
      char name[5], res[4];
      strncpy(name, line + 12, 4); name[4] = '\0';
      strncpy(res, line + 17, 3); res[3] = '\0';
      
      // Trim whitespace
      char *clean_name = strtok(name, " ");
      
      if (is_sidechain(res, clean_name)) {
        atoms[count].atom_no = atoi(line + 6);
        strcpy(atoms[count].atom_name, clean_name);
        strcpy(atoms[count].res_name, res);
        atoms[count].res_no = atoi(line + 22);
        atoms[count].x = atof(line + 30);
        atoms[count].y = atof(line + 38);
        atoms[count].z = atof(line + 46);
        count++;
      }
    }
  }
  fclose(infile);
  
  // 2. Global Minimum Variables
  int i,j;
  double min_dist = 1e9; 
  int atom1_idx = -1, atom2_idx = -1;
  
  // 3. Matrix Header likhna (CSV ke liye)
  fprintf(outfile, "Atom_ID");
  for (i = 0; i < count; i++) {
    fprintf(outfile, ",%s_%d_%s", atoms[i].atom_name, atoms[i].res_no, atoms[i].res_name);
  }
  fprintf(outfile, "\n");
  
  // 4. Nested Loops for Distance Calculation
  for (i = 0; i < count; i++) {
    fprintf(outfile, "%s_%d_%s", atoms[i].atom_name, atoms[i].res_no, atoms[i].res_name);
    
    for (j = 0; j < count; j++) {
      if (atoms[i].res_no == atoms[j].res_no) {
        // Same residue ke beech distance 0
        fprintf(outfile, ",0.000");
      } else {
        double dx = atoms[i].x - atoms[j].x;
        double dy = atoms[i].y - atoms[j].y;
        double dz = atoms[i].z - atoms[j].z;
        double dist = sqrt(dx*dx + dy*dy + dz*dz);
        
        fprintf(outfile, ",%.3f", dist);
        
        // Global Minimum Update
        if (dist >2.5){
          if (dist < min_dist) {
            min_dist = dist;
            atom1_idx = i;
            atom2_idx = j;
          }
        }
      }
    }
    fprintf(outfile, "\n");
  }
  
  fclose(outfile);
  
  // 5. Result Display
  if (atom1_idx != -1) {
    printf("\n--- Global Minimum Distance Result ---\n");
    printf("Min Distance: %.4f Angstroms\n", min_dist);
    printf("Atom 1: %s (No: %d) in Residue %s (No: %d)\n", 
           atoms[atom1_idx].atom_name, atoms[atom1_idx].atom_no, 
           atoms[atom1_idx].res_name, atoms[atom1_idx].res_no);
    printf("Atom 2: %s (No: %d) in Residue %s (No: %d)\n", 
           atoms[atom2_idx].atom_name, atoms[atom2_idx].atom_no, 
           atoms[atom2_idx].res_name, atoms[atom2_idx].res_no);
    printf("\nExcel file 'sidechain_matrix.csv' generate ho gayi hai.\n");
  }
  free(atoms);
  return 0;
}