#include "report.h"
#include "stdio.h"

void record_report (const char *file)
{
    FILE *fp;
    if (!(fp = fopen (file, "w")))
    {
        printf ("Cannot open file for write report!\n");
        return;
    }

    // table with result

    fclose (fp);
}

void record_file (const char *file, double *H, int mh, double *V, int mv, int mx)
{
    FILE *fp;
    int i = 0;

    if (!(fp = fopen (file, "w")))
    {
        printf ("Cannot open file for write data!\n");
        return;
    }

    fprintf (fp, "\n\n\n>>>>>> Density >>>>>>>\n Size = %d\n", mh);
    for (i = 0; i < mh; i++)
    {
        fprintf (fp, " %f ", H[i]);
        if (i % (mx - 1) == mx - 2)
            fprintf (fp, "\n");
    }
    fprintf (fp, "\n\n\n<<<<<< Density <<<<<<<");

    fprintf (fp, "\n\n\n>>>>>> Pulse >>>>>>>\n Size = %d\n", mv);
    fprintf (fp, "I\n");
    for (i = 0; i < mh; i++)
    {
        fprintf (fp, " %f ", V[i]);
        if (i % mx == mx - 1)
            fprintf (fp, "\n");
    }
    fprintf (fp, "II\n");
    for (; i < mh; i++)
    {
        fprintf (fp, " %f ", V[i]);
        if (i % mx == mx - 1)
            fprintf (fp, "\n");
    }

    fprintf (fp, "\n\n\n<<<<<< Pulse <<<<<<<");

    fclose (fp);
}
