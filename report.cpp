#include "report.h"

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