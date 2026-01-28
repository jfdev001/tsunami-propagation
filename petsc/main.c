static char help[] = "Hello world for PETSC\n\n";

#include "stdio.h"
#include <petscvec.h>

int main(int argc, char **args) {
  printf("Hello world!\n");
  PetscFunctionBeginUser;
  PetscCall(PetscInitialize(&argc, &args, NULL, help));
  PetscCall(PetscFinalize());
  return 0;
}
