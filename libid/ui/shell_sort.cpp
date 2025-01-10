// SPDX-License-Identifier: GPL-3.0-only
//
#include "ui/shell_sort.h"

#include "misc/prototyp.h" // for stricmp

void shell_sort(void *v1, int n, unsigned sz)
{
    auto lc_compare = [](void *arg1, void *arg2) // for sort
    {
        char **choice1 = (char **) arg1;
        char **choice2 = (char **) arg2;

        return stricmp(*choice1, *choice2);
    };

    char *v = (char *) v1;
    for (int gap = n/2; gap > 0; gap /= 2)
    {
        for (int i = gap; i < n; i++)
        {
            for (int j = i-gap; j >= 0; j -= gap)
            {
                if (lc_compare((char **)(v+j*sz), (char **)(v+(j+gap)*sz)) <= 0)
                {
                    break;
                }
                void *temp = *(char **) (v + j * sz);
                *(char **)(v+j*sz) = *(char **)(v+(j+gap)*sz);
                *(char **)(v+(j+gap)*sz) = static_cast<char *>(temp);
            }
        }
    }
}
