"""
Example usage of xyzPlots module
"""
import xyzPlots
import pandas as pd


def main():
    print('Reading in data...')
    df = pd.read_csv('data/2hdme_tiny_data.csv')

    x = df.Lambda2
    y = df.Lambda4
    z = df.uv_ef
    xyzPlots.Heatmap(x, y, z, xLabel=r'$\lambda_2$', yLabel=r'$\lambda_4$',
                     zLabel=r'max $\Lambda$', output='plots/example.png')

    print('Plot saved in plots/example.pdf')


if __name__ == "__main__":
    main()
