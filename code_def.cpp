#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <fstream>
#include <sstream>
using namespace std;

int main()
{
  bool x;
  string vib, line, bec, mass, check;
  int cellsize, n = 0, I = 0, j = 0;
  cout << "Input cell size \n";
  cin >> cellsize;
  cout << "Input vibrations file (.txt) \n";
  cin >> vib;
  cout << "Input masses file (.txt) \n";
  cin >> mass;
  cout << "Input Born Effective Charge file (.txt) \n";
  cin >> bec;
  cout << "Input reference file for Harmonic Oscillator Strenght (.txt) \n";
  cin >> check;

  double vibrations[3 * cellsize][3 * cellsize];
  double masses[cellsize];
  string labels[cellsize];
  double borneff[cellsize][3][3];
  double oscstr[3][3 * cellsize];
  double ref[3 * cellsize];
  double delta[3 * cellsize];
  double HOS [3 * cellsize];

  ifstream myfile(vib);
  if (myfile.is_open())
  {
    while (getline(myfile, line))
    {
      stringstream iss(line);
      double a, b, c, d, e, f;
      if (!(iss >> a >> b >> c >> d >> e >> f))
      {
        break;
      }
      vibrations[I][j] = a;
      vibrations[I + 1][j] = c;
      vibrations[I + 2][j] = e;
      I = I + 3;

      if (I == 3 * cellsize)
      {
        I = 0;
        j++;
      }
    }
    myfile.close();
  }
  else
  {
    cout << "Unable to open file " << vib;
    return 1;
  }

  cout << "Do you want to check the file " << vib << " ? [Y=1;N=0] \n";
  cin >> x;
  if (x == 1)
  {
    cout << "CHECK INPUT FOR FILE " << vib << '\n';
    for (int i = 0; i < 3 * cellsize; i++)
    {
      cout << "vibrational mode " << i + 1 << '\n';
      for (int n = 0; n < cellsize; n++)
      {
        cout << vibrations[3 * n][i] << '\t' << vibrations[3 * n + 1][i] << '\t' << vibrations[3 * n + 2][i] << '\n';
      }
    }
  }

  j = 0;

  ifstream ourfile(mass);
  if (ourfile.is_open())
  {
    while (getline(ourfile, line))
    {
      stringstream iss(line);
      double a, c;
      string b;
      if (!(iss >> a >> b >> c))
      {
        break;
      }
      labels[j] = b;
      masses[j] = c;
      j++;
    }
  }
  else
  {
    cout << "Unable to open file " << mass;
    return 1;
  }

  cout << "Do you want to check the file " << mass << " ? [Y=1;N=0] \n";
  cin >> x;
  if (x == 1)
  {
    cout << "CHECK INPUT FOR FILE " << mass << '\n';
    for (int n = 0; n < cellsize; n++)
    {
      cout << labels[n] << '\t' << masses[n] << '\t' << sqrt(masses[n]) << '\n';
    }
    cout << '\n';
  }

  j = 0;
  I = 0;

  ifstream themfile(bec);
  if (themfile.is_open())
  {

    while (getline(themfile, line))
    {
      stringstream iss(line);
      double a, b, c;
      if (!(iss >> a >> b >> c))
      {
        break;
      }
      borneff[I][j][0] = a;
      borneff[I][j][1] = b;
      borneff[I][j][2] = c;
      j++;
      if (j == 3)
      {
        j = 0;
        I++;
      }
    }
  }
  else
  {
    cout << "Unable to open file " << bec;
    return 1;
  }

  for (int i = 0; i < 3; i++)
  {
    for (int j = 0; j < 3; j++)
    {
      double sum = 0;
      for (int n = 0; n < cellsize; n++)
      {
        sum += borneff[n][i][j];
      }
      for (int n = 0; n < cellsize; n++)
      {
        borneff[n][i][j] -= (sum / cellsize);
      }
    }
  }

  cout << "BEC tensors corrected for acoustic sum rule \n";

  cout << "Do you want to check the file " << bec << " ? [Y=1;N=0] \n";
  cin >> x;
  if (x == 1)
  {
    cout << "CHECK INPUT FOR FILE " << vib << '\n';
    for (int i = 0; i < cellsize; i++)
    {
      cout << "BEC tensor for atom " << i + 1 << '\n';
      for (int n = 0; n < 3; n++)
      {
        cout << borneff[i][n][0] << '\t' << borneff[i][n][1] << '\t' << borneff[i][n][2] << '\n';
      }
    }
  }

  ifstream hefile(check);
  if (hefile.is_open())
  {
    while (getline(hefile, line))
    {
      stringstream iss(line);
      double a, b, c, d;
      if (!(iss >> a >> b >> c >> d))
      {
        break;
      }
      ref[j] = d;
      j++;
    }
  }
  else
  {
    cout << "Unable to open file " << check;
    return 1;
  }

  cout << "Correcting and renormalising vibration vectors \n";
  cout << "Do you want to check the renormalisation? [Y=1;N=0] \n";
  cin >> x;

  for (int i = 0; i < 3 * cellsize; i++)
  {

    double lenght = 0;
    for (int n = 0; n < cellsize; n++)
    {

      double k = vibrations[3 * n][i] * vibrations[3 * n][i] +
                 vibrations[3 * n + 1][i] * vibrations[3 * n + 1][i] +
                 vibrations[3 * n + 2][i] * vibrations[3 * n + 2][i];

      lenght += k;
    }
    double norm = sqrt(lenght);

    if (x == 1)
    {
      cout << "The norm of v_In as input for the mode " << i + 1 << " is " << norm << '\n';
    }

    for (int n = 0; n < cellsize; n++)
    {
      // cout << sqrt(masses[n]) << '\n';
      vibrations[3 * n][i] = vibrations[3 * n][i] * sqrt(masses[n]);
      vibrations[3 * n + 1][i] = vibrations[3 * n + 1][i] * sqrt(masses[n]);
      vibrations[3 * n + 2][i] = vibrations[3 * n + 2][i] * sqrt(masses[n]);
    }

    double lenght2 = 0;
    for (int n = 0; n < cellsize; n++)
    { // indice n scorre gli atomi

      double k = vibrations[3 * n][i] * vibrations[3 * n][i] +
                 vibrations[3 * n + 1][i] * vibrations[3 * n + 1][i] +
                 vibrations[3 * n + 2][i] * vibrations[3 * n + 2][i];

      lenght2 += k;
    }
    double norm2 = sqrt(lenght2);

    if (x == 1)
    {
      cout << "The norm of the vector u_In= v_In * sqrt(M_I) for the mode" << i + 1 << " is " << norm2 << '\n';
    }

    for (int n = 0; n < cellsize; n++)
    {
      vibrations[3 * n][i] /= norm2;
      vibrations[3 * n + 1][i] /= norm2;
      vibrations[3 * n + 2][i] /= norm2;
    }

    double lenght3 = 0;
    for (int n = 0; n < cellsize; n++)
    { // indice n scorre gli atomi

      double k = vibrations[3 * n][i] * vibrations[3 * n][i] +
                 vibrations[3 * n + 1][i] * vibrations[3 * n + 1][i] +
                 vibrations[3 * n + 2][i] * vibrations[3 * n + 2][i];

      lenght3 += k;
    }

    double norm3 = sqrt(lenght3);

    if (x == 1)
    {
      cout << "The renormalised u_In norm for the mode" << i + 1 << " is " << norm3 << '\n';
    }
  }

  cout << "Do you want to ckeck orthogonality for u_n? [Y=1;N=0] \n";
  cin >> x;

  for (int i = 0; i < 3 * cellsize; i++)
  {

    int j = 3 * cellsize - 1;
    while (j > i)
    {
      double ort = 0;
      for (int n = 0; n < cellsize; n++)
      {
        double diff = vibrations[3 * n][i] * vibrations[3 * n][j] +
                      vibrations[3 * n + 1][i] * vibrations[3 * n + 1][j] +
                      vibrations[3 * n + 2][i] * vibrations[3 * n + 2][j];
        ort += diff;
      }

      if (x == 1)
      {
        cout << "Mode " << i + 1 << " scalar " << j + 1 << " = " << ort << '\n';
      }

      j = j - 1;
    }
  }

  cout << "Computing Harmonic Oscillator Strenghts \n";

  for (int i = 0; i < 3 * cellsize; i++)
  {
    for (int k = 0; k < 3; k++)
    {
      double s = 0;
      for (int n = 0; n < cellsize; n++)
      {
        for (int j = 0; j < 3; j++)
        {
          s += (vibrations[3 * n + j][i] * borneff[n][k][j]) / sqrt(masses[n]);
        }
      }
      oscstr[k][i] = s;
    }
    HOS[i] = pow(oscstr[0][i], 2) + pow(oscstr[1][i], 2) + pow(oscstr[2][i], 2);
    cout << "HOS for mode " << i + 1 << " is " << HOS[i] << '\n';
  }

  cout << "Do you want to calculate the corrective factor or do you want to manually input it? [1/0] \n";
  cin >> x;
  double mean = 0;
  if (x == 1)
  {
    double s = 0;
    for (int i = 0; i < 3 * cellsize; i++)
    {
      double corr = ref[i] / HOS[i];
      s += corr;
      cout << "The corrective factor for the mode " << i + 1 << " is " << corr << '\n';
    }
    mean = s / (3 * cellsize - 3);
    cout << "The average corrective factor is " << mean << '\n';
  }

  if (x == 0)
  {
    cout << "Input corrective factor \n";
    cin >> mean;
  }

  cout << "Do you want to disable vibrations for certain atoms? [Y=1;N=0] \n";
  bool y;
  cin >> x;
  if (x == 1)
  {

    cout << "Do you want to chose by atom name or by atom index inside the crystal? [1=name;0=index] \n";
    cin >> y;
    if (y == 1)
    {
      string name;
      int I = 0;
      while (I < 1)
      {
        cout << "Input atom name, type 'out' to exit the selection \n";
        cin >> name;
        if (name == "out")
        {
          I++;
        }
        for (int n = 0; n < cellsize; n++)
        {
          if (labels[n] == name)
          {
            for (int i = 0; i < 3 * cellsize; i++)
            {
              vibrations[3 * n][i] = 0.;
              vibrations[3 * n + 1][i] = 0.;
              vibrations[3 * n + 2][i] = 0.;
            }
            cout << "Deleted vibrations for the atom " << labels[n] << " with index " << n + 1 << '\n';
          }
        }
      }
      cout << "Do you want to check the edited vibrations file? [Y=1;N=0] \n";
      cin >> x;
      if (x == 1)
      {
        cout << "CHECK INPUT FOR FILE " << vib << '\n';
        for (int i = 0; i < 3 * cellsize; i++)
        {
          cout << "vibrational mode " << i + 1 << '\n';
          for (int n = 0; n < cellsize; n++)
          {
            cout << vibrations[3 * n][i] << '\t' << vibrations[3 * n + 1][i] << '\t' << vibrations[3 * n + 2][i] << '\n';
          }
        }
      }
    }

    if (y == 0)
    {
      int n;
      int I = 0;
      while (I < 1)
      {
        cout << "Input atom index [1;" << cellsize << "], type '0' to exit the selection \n";
        cin >> n;
        bool y = 1;
        if (n == 0)
        {
          y = 0;
          I++;
        }
        if (y == 1)
        {
          for (int i = 0; i < 3 * cellsize; i++)
          {
            vibrations[3 * (n - 1)][i] = 0.;
            vibrations[3 * (n - 1) + 1][i] = 0.;
            vibrations[3 * (n - 1) + 2][i] = 0.;
          }
          cout << "Deleted vibrations for the atom " << labels[n - 1] << " with index " << n << '\n';
        }
      }
      cout << "Do you want to check the edited vibrations file? [Y=1;N=0] \n";
      cin >> x;
      if (x == 1)
      {
        cout << "CHECK INPUT FOR FILE " << vib << '\n';
        for (int i = 0; i < 3 * cellsize; i++)
        {
          cout << "vibrational mode " << i + 1 << '\n';
          for (int n = 0; n < cellsize; n++)
          {
            cout << vibrations[3 * n][i] << '\t' << vibrations[3 * n + 1][i] << '\t' << vibrations[3 * n + 2][i] << '\n';
          }
        }
      }
    }
  }

  cout << "Computing delta correction \n";
  for (int i = 0; i < 3 * cellsize; i++)
  {
    double norm = 0;
    for (int n = 0; n < cellsize; n++)
    {
      norm += vibrations[3 * n][i] * vibrations[3 * n][i] +
              vibrations[3 * n + 1][i] * vibrations[3 * n + 1][i] +
              vibrations[3 * n + 2][i] * vibrations[3 * n + 2][i];
    }
    delta[i] = norm;
  }
  cout << "Do you want to check the delta corrections? [Y=1;N=0] \n";
  cin >> x;
  if (x == 1)
  {
    for (int i = 0; i < 3 * cellsize; i++)
    {
      cout << "vibrational mode " << i + 1 << '\t' << delta[i] << '\n';
    }
  }

  cout << "Corrected HOS (with delta and corrective factor) are \n";

  for (int i = 0; i < 3 * cellsize; i++)
  {
    HOS[i] = HOS[i] * delta[i] * mean;
    cout << "Mode " << i + 1 << '\t'<< HOS[i] << '\n';
  }

  return 0;
}