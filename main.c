
/*********************************************************************************************

    This is public domain software that was developed by or for the U.S. Naval Oceanographic
    Office and/or the U.S. Army Corps of Engineers.

    This is a work of the U.S. Government. In accordance with 17 USC 105, copyright protection
    is not available for any work of the U.S. Government.

    Neither the United States Government, nor any employees of the United States Government,
    nor the author, makes any warranty, express or implied, without even the implied warranty
    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE, or assumes any liability or
    responsibility for the accuracy, completeness, or usefulness of any information,
    apparatus, product, or process disclosed, or represents that its use would not infringe
    privately-owned rights. Reference herein to any specific commercial products, process,
    or service by trade name, trademark, manufacturer, or otherwise, does not necessarily
    constitute or imply its endorsement, recommendation, or favoring by the United States
    Government. The views and opinions of authors expressed herein do not necessarily state
    or reflect those of the United States Government, and shall not be used for advertising
    or product endorsement purposes.

*********************************************************************************************/

 /********************************************************************
 *
 * Module Name : main.c
 *
 * Author/Date : Jan C. Depner - 10/21/16
 *
 * Description : Filters out returns that have digitizer noise or 
 *               whose starting waveform amplitude exceeds the user
 *               defined threshold
 *
 ********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <errno.h>
#include <string.h>
#include <getopt.h>
#include <math.h>


/* Local Includes. */

#include "nvutility.h"

#include "czmil.h"

#include "version.h"


void usage ()
{
  fprintf (stderr, "\nUsage: czmil_noise_filter [-1] [-2] [-3] [-4] [-5] [-6] [-7] [-9] [-a THRESHOLD] [-s SHAL_AMP] [-d DEEP_AMP] CZMIL_CPF_FILENAME\n");
  fprintf (stderr, "Where:\n");
  fprintf (stderr, "\t-1 = filter channel 1\n");
  fprintf (stderr, "\t-2 = filter channel 2\n");
  fprintf (stderr, "\t-3 = filter channel 3\n");
  fprintf (stderr, "\t-4 = filter channel 4\n");
  fprintf (stderr, "\t-5 = filter channel 5\n");
  fprintf (stderr, "\t-6 = filter channel 6\n");
  fprintf (stderr, "\t-7 = filter channel 7\n");
  fprintf (stderr, "\t-9 = filter channel 9\n");
  fprintf (stderr, "\tTHRESHOLD = waveform amplitude second difference change threshold [default = noise filter disabled]\n");
  fprintf (stderr, "\tSHAL_AMP = Shallow channel starting amplitude threshold [default = shallow amplitude filter disabled]\n");
  fprintf (stderr, "\tDEEP_AMP = Deep channel starting amplitude threshold [default = deep amplitude filter disabled]\n\n");
  exit (-1);
}



int32_t main (int32_t argc, char **argv)
{
  char               cpf_file[1024], cwf_file[1024];
  int32_t            cpf_hnd = -1, cwf_hnd = -1, i, j, k, m, threshold = 0, length, percent = 0, old_percent = -1, kill_count = 0, shal_amp = 0, deep_amp = 0, check_filt;
  int16_t            *diff;
  CZMIL_CPF_Header   cpf_header;
  CZMIL_CPF_Data     cpf;
  CZMIL_CWF_Header   cwf_header;
  CZMIL_CWF_Data     cwf;
  uint8_t            channel[9] = {NVFalse, NVFalse, NVFalse, NVFalse, NVFalse, NVFalse, NVFalse, NVFalse, NVFalse}, check_it, mod_flag;
  char               c;
  extern char        *optarg;
  extern int         optind;


  fprintf (stderr, "\n\n %s \n\n\n", VERSION);


  while ((c = getopt (argc, argv, "12345679a:s:d:")) != EOF)
    {
      switch (c)
        {
        case 'a':
          sscanf (optarg, "%d", &threshold);
          break;

        case 's':
          sscanf (optarg, "%d", &shal_amp);
          break;

        case 'd':
          sscanf (optarg, "%d", &deep_amp);
          break;

        case '1':
          channel[0] = NVTrue;
          break;

        case '2':
          channel[1] = NVTrue;
          break;

        case '3':
          channel[2] = NVTrue;
          break;

        case '4':
          channel[3] = NVTrue;
          break;

        case '5':
          channel[4] = NVTrue;
          break;

        case '6':
          channel[5] = NVTrue;
          break;

        case '7':
          channel[6] = NVTrue;
          break;

        case '9':
          channel[8] = NVTrue;
          break;

        default:
          usage ();
          break;
        }
    }


  /*  Make sure we got the mandatory file name argument, at least one channel, and one kind of filter.  */

  check_it = NVFalse;
  for (i = 0 ; i < 9 ; i++)
    {
      if (channel[i])
        {
          check_it = NVTrue;
          break;
        }
    }

  check_filt = threshold + shal_amp + deep_amp;

  if (!check_it || check_filt <= 0 || optind >= argc) usage ();


  strcpy (cpf_file, argv[optind]);


  if (!strstr (cpf_file, ".cpf")) usage ();


  if ((cpf_hnd = czmil_open_cpf_file (cpf_file, &cpf_header, CZMIL_UPDATE)) < 0)
    {
      czmil_perror ();
      exit (-1);
    }


  strcpy (cwf_file, cpf_file);
  sprintf (&cwf_file[strlen (cwf_file) - 4], ".cwf");

  if ((cwf_hnd = czmil_open_cwf_file (cwf_file, &cwf_header, CZMIL_READONLY)) < 0)
    {
      czmil_perror ();
      exit (-1);
    }


  fprintf (stderr, "\n\n File : %s\n\n", cpf_file);


  for (i = 0 ; i < cpf_header.number_of_records ; i++)
    {
      if (czmil_read_cpf_record (cpf_hnd, i, &cpf) != CZMIL_SUCCESS)
        {
          czmil_perror ();
          exit (-1);
        }


      if (czmil_read_cwf_record (cwf_hnd, i, &cwf) != CZMIL_SUCCESS)
        {
          czmil_perror ();
          exit (-1);
        }


      mod_flag = NVFalse;

      for (j = 0 ; j < 8 ; j++)
        {
          if (channel[j])
            {
              check_it = NVFalse;


              /*  We need to check for valid returns.  */

              for (k = 0 ; k < cpf.returns[j] ; k++)
                {
                  /*  If we are doing a digitizer noise filter, we need to change old CZMIL_RETURN_FILTER_INVAL/CZMIL_DIGITIZER_NOISE points
                      back to valid.  */

                  if (threshold > 0)
                    {
                      if (cpf.channel[j][k].status & CZMIL_RETURN_FILTER_INVAL)
                        {
                          if (cpf.channel[j][k].filter_reason == CZMIL_DIGITIZER_NOISE)
                            {
                              cpf.channel[j][k].status &= ~CZMIL_RETURN_FILTER_INVAL;
                              cpf.channel[j][k].filter_reason = CZMIL_WAVEFORM_VALID;
                              mod_flag = NVTrue;
                            }
                        }
                    }

                  if (!(cpf.channel[j][k].status & CZMIL_RETURN_INVAL)) check_it = NVTrue;
                }


              /*  No point in checking if there are no valid returns.  */

              if (check_it)
                {
                  /*  Check for digitizer noise if requested.  */

                  if (threshold > 0)
                    {
                      length = cwf.number_of_packets[j] * 64;


                      diff = (int16_t *) calloc (length - 1, sizeof (int16_t));
                      if (diff == NULL)
                        {
                          perror ("Allocating diff memory");
                          exit (-1);
                        }


                      for (k = 1 ; k < length ; k++) diff[k - 1] = cwf.channel[j][k] - cwf.channel[j][k - 1];


                      for (k = 1 ; k < length - 1 ; k++)
                        {
                          if (diff[k] - diff[k - 1] > threshold)
                            {
                              for (m = 0 ; m < cpf.returns[j] ; m++)
                                {
                                  if (!(cpf.channel[j][m].status & CZMIL_RETURN_INVAL))
                                    {
                                      cpf.channel[j][m].status |= CZMIL_RETURN_FILTER_INVAL;
                                      cpf.channel[j][m].filter_reason = CZMIL_DIGITIZER_NOISE;
                                      kill_count++;
                                      mod_flag = NVTrue;
                                    }
                                }
                            }
                        }

                      free (diff);
                    }


                  /*  Check the channel for a starting amplitude higher than the amplitude threshold (if requested)  */

                  check_filt = 0;
                  if (j == 8)
                    {
                      if (deep_amp > 0)
                        {
                          if (cwf.channel[j][0] > deep_amp) check_filt = 1;
                        }
                    }
                  else
                    {
                      if (shal_amp > 0)
                        {
                          if (cwf.channel[j][0] > shal_amp) check_filt = 1;
                        }
                    }

                  if (check_filt)
                    {
                      for (m = 0 ; m < cpf.returns[j] ; m++)
                        {
                          if (!(cpf.channel[j][m].status & CZMIL_RETURN_INVAL))
                            {
                              cpf.channel[j][m].status |= CZMIL_RETURN_FILTER_INVAL;
                              cpf.channel[j][m].filter_reason = CZMIL_START_AMP_EXCEEDS_THRESHOLD;
                              kill_count++;
                              mod_flag = NVTrue;
                            }
                        }
                    }
                }
            }
        }


      if (mod_flag)
        {
          if (czmil_update_cpf_return_status (cpf_hnd, i, &cpf) != CZMIL_SUCCESS)
            {
              czmil_perror ();
              exit (-1);
            }
        }


      percent = NINT (((float) i / (float) cpf_header.number_of_records) * 100.0);
      if (old_percent != percent)
        {
          fprintf (stdout, "%3d%% processed    \r", percent);
          fflush (stdout);
          old_percent = percent;
        }
    }


  fprintf (stdout, "100%% processed, %d invalidated\n", kill_count);
  fflush (stdout);

  czmil_close_cwf_file (cwf_hnd);
  czmil_close_cpf_file (cpf_hnd);

  return (0);
}
