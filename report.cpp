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

void record_file(const char *file)
{
    FILE *fp;
    if (!(fp = fopen (file, "w")))
    {
        printf ("Cannot open file for write report!\n");
        return;
    }

    // data

    fclose (fp);
}
