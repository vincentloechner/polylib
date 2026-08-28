#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <polylib/polylib.h>

/****************************************************/
/* Read_ParamNames() :                              */
/* Reads FILE *in for the parameter names           */
/* if in==NULL or not enough parameters on input,   */
/*  use default names                               */
/* returns an n-array of strings                    */
/****************************************************/
const char **Read_ParamNames(FILE *in, int m) {
  int c, i, j, f;
  char s[1024], param[32];

  if (!in)
    f = 0;
  else
    do
      f = (fgets(s, 1024, in) != NULL);
    while (f && (*s == '#' || *s == '\n'));

  const char **param_name = malloc(m * sizeof(char *));
  i = 0;
  if (f) {
    c = 0;
    for (; i < m; ++i) {
      j = 0;
      for (;; ++c) {
        if (s[c] == ' ') {
          if (j == 0)
            continue;
          else
            break;
        }
        if (s[c] == '\n' || s[c] == 0)
          break;
        param[j++] = s[c];
      }

      /* Not enough parameters on input, end */
      if (j == 0)
        break;
      param[j] = 0;
      char *name = malloc((j + 1) * sizeof(char));
      strcpy(name, param);
      param_name[i] = name;
    }
  }

  /* Not enough parameters on input : use default names */
  if (!f || i != m) {
    for (; i < m; ++i) {
      char *name = malloc(2 * sizeof(char));
      snprintf(name, 2, "%c", PCHAR + i + 1);
      param_name[i] = name;
    }
  }
  return param_name;
} /* Read_ParamNames */

void Free_ParamNames(const char **params, int m) {
  while (--m >= 0) {
    // only exception modify the constant pointer, by freeing it at the end...
    free((void *)params[m]);
  }
  free(params);
}
