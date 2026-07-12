import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde


def _pairwise_complete_obs(x, y):
    """Return x, y with pairwise complete observations (no NaNs in either)."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    mask = ~np.isnan(x) & ~np.isnan(y)
    return x[mask], y[mask]


def calc_quadrant_analysis(
    xval,
    yval,
    do_normalization=True,
    hole_sizes=range(0, 11),
    orient="+",
    plot=True,
    **plot_kwargs,
):
    """
    Calculating Coherent Structures following Quadrant Analysis (Python version)

    Parameters
    ----------
    xval : array-like
        Values of x variable.
    yval : array-like
        Values of y variable.
    do_normalization : bool, default True
        Normalize values: (x - mean(x))/sd(x), ignoring NaNs.
    hole_sizes : sequence of int, default range(0, 11)
        Desired hole sizes (integers >= 0).
    orient : {"+", "-"}, default "+"
        Orientation for down-gradient flux.
    plot : bool, default True
        If True, calls plot_quadrant_analysis.
    **plot_kwargs
        Passed to plot_quadrant_analysis.

    Returns
    -------
    dict with keys:
        "hole_sizes", "occurrence", "product", "covariance",
        "covariance_total", "correlation_total",
        "exuberance", "organization_ratio", "meta"
    """

    xval = np.asarray(xval, dtype=float)
    yval = np.asarray(yval, dtype=float)
    hole_sizes = np.asarray(list(hole_sizes), dtype=float)
    x_clean, y_clean = _pairwise_complete_obs(xval, yval)
    if x_clean.size > 1:
        covariance_total = np.cov(x_clean, y_clean)[0, 1]
        correlation_total = np.corrcoef(x_clean, y_clean)[0, 1]
    else:
        covariance_total = np.nan
        correlation_total = np.nan
    if do_normalization:
        x_mean = np.nanmean(xval)
        x_std = np.nanstd(xval,ddof=1)
        y_mean = np.nanmean(yval)
        y_std = np.nanstd(yval,ddof=1)
        xval = (xval-x_mean)/x_std
        yval = (yval-y_mean)/y_std
    nh = len(hole_sizes)
    occurrence = np.full((4,nh),np.nan)
    product = np.full((4,nh),np.nan)
    covariance = np.full((4,nh),np.nan)
    prod = xval*yval
    abs_prod = np.abs(prod)
    #holes
    for i, hole in enumerate(hole_sizes):
        is_q1 = (xval>0) & (yval>0) & (abs_prod>hole)
        is_q2 = (xval<0) & (yval>0) & (abs_prod>hole)
        is_q3 = (xval<0) & (yval<0) & (abs_prod>hole)
        is_q4 = (xval>0) & (yval<0) & (abs_prod>hole)
        masks = [is_q1, is_q2, is_q3, is_q4]
        for q in range(4):
            m = masks[q]
            if np.any(m):
                xv_q = xval[m]
                yv_q = yval[m]
                occurrence[q,i] = np.sum(m)
                product[q,i] = np.nanmean(xv_q*yv_q)
                xv_q_clean, yv_q_clean = _pairwise_complete_obs(xv_q,yv_q)
                if xv_q_clean.size > 1:
                    covariance[q,i] = np.cov(xv_q_clean,yv_q_clean)[0,1]
                else:
                    covariance[q,i] = np.nan
            else:
                occurrence[q,i] = 0
                product[q,i] = np.nan
                covariance[q,i] = np.nan
    #diagnostics: exuberance and organization ratio
    if orient == "-":
        exub = (product[0,:]+product[2,:])/(product[1,:]+product[3,:])
        oratio = (occurrence[0,:]+occurrence[2,:])/(occurrence[1,:]+occurrence[3,:])
    elif orient == "+":
        exub = (product[1, :]+product[3,:])/(product[0,:]+product[2,:])
        oratio = (occurrence[1,:]+occurrence[3,:])/(occurrence[0,:]+occurrence[2,:])
    else:
        raise ValueError('The orientation has to be either "+" or "-".')

    qa_out = {
        "hole_sizes": hole_sizes,
        "occurrence": occurrence,
        "product": product,
        "covariance": covariance,
        "covariance_total": covariance_total,
        "correlation_total": correlation_total,
        "exuberance": exub,
        "organization_ratio": oratio,
        "meta": (
            "Output format: rows represent the quadrants Q1, Q2, Q3, Q4 -- "
            "columns represent selected hole sizes"
        ),
    }
    if plot:
        plot_quadrant_analysis(xval, yval, do_normalization=False, **plot_kwargs)
        print(qa_out)
    return qa_out


def plot_quadrant_analysis(
    xval,
    yval,
    do_normalization=True,
    hole_sizes=(1, 2, 3),
    plot_kde2d=True,
    contours=None,
    print_fit=True,
    **plot_kwargs,
):
    """
    Plotting Quadrant Analysis (Python version)

    Parameters
    ----------
    xval, yval : array-like
        Values of variables to plot.
    do_normalization : bool, default True
        Normalize values: (x - mean(x))/sd(x).
    hole_sizes : sequence of int, default (1, 2, 3)
        Hole sizes for hyperbolic curves.
    plot_kde2d : bool, default True
        If True, plot 2D KDE contour lines.
    contours : sequence of float, optional
        Levels of contour lines; default 10**(-3:3) as in R.
    print_fit : bool, default True
        If True, print linear regression summary.
    **plot_kwargs
        Passed to plt.scatter.
    """

    xval = np.asarray(xval, dtype=float)
    yval = np.asarray(yval, dtype=float)
    hole_sizes = np.asarray(list(hole_sizes), dtype=float)
    if do_normalization:
        x_mean = np.nanmean(xval)
        x_std = np.nanstd(xval,ddof=1)
        y_mean = np.nanmean(yval)
        y_std = np.nanstd(yval,ddof=1)
        xval = (xval-x_mean)/x_std
        yval = (yval-y_mean)/y_std
    #scatter plot
    col = (0.6,0.6,0.6,0.1)
    plt.scatter(xval,yval,c=[col],s=10,**plot_kwargs)
    plt.axhline(0,linestyle="--",color="black",linewidth=2)
    plt.axvline(0,linestyle="--",color="black",linewidth=2)
    mask = ~np.isnan(xval) & ~np.isnan(yval)
    x_clean = xval[mask]
    y_clean = yval[mask]
    if x_clean.size > 1:
        X = np.vstack([np.ones_like(x_clean),x_clean]).T
        beta, _, _, _ = np.linalg.lstsq(X,y_clean,rcond=None)
        a,b = beta
        if print_fit:
            y_pred = a+b*x_clean
            residuals = y_clean-y_pred
            ss_res = np.sum(residuals**2)
            ss_tot = np.sum((y_clean-np.mean(y_clean))**2)
            r2 = 1-ss_res/ss_tot if ss_tot>0 else np.nan
            print("Linear regression fit:")
            print(f"  intercept (a): {a:.4f}")
            print(f"  slope (b):     {b:.4f}")
            print(f"  R^2:           {r2:.4f}")
        x_line = np.linspace(np.nanmin(xval),np.nanmax(xval),200)
        y_line = a+b*x_line
        plt.plot(x_line,y_line, color="darkred", linewidth=3)
    else:
        if print_fit:
            print("Not enough non-NaN points for linear regression.")
    if plot_kde2d and x_clean.size > 1:
        if contours is None:
            contours = 10.0**np.arange(-3,4)  # 10^(-3:3)
        contours = np.asarray(contours,dtype=float)
        nc = len(contours)
        cmap = plt.get_cmap("coolwarm")
        colors = [cmap(i/(nc-1)) for i in range(nc)]
        kde = gaussian_kde(np.vstack([x_clean,y_clean]))
        xmin, xmax = np.nanmin(x_clean),np.nanmax(x_clean)
        ymin, ymax = np.nanmin(y_clean),np.nanmax(y_clean)
        xx, yy = np.meshgrid(
            np.linspace(xmin,xmax,100),
            np.linspace(ymin,ymax,100),
        )
        positions = np.vstack([xx.ravel(),yy.ravel()])
        zz = kde(positions).reshape(xx.shape)
        plt.contour(xx,yy,zz,levels=contours,colors=colors,linewidths=1)
    xs = np.linspace(-10,10,400)
    xs = xs[xs!=0.0]
    for h in hole_sizes:
        ys_pos = h/xs
        ys_neg = -h/xs
        plt.plot(xs,ys_pos,color="black",linestyle="--",linewidth=1)
        plt.plot(xs,ys_neg,color="black",linestyle="--",linewidth=1)
    plt.legend(["linear fit","hyperbolic holes"],loc="lower left")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.title("Quadrant Analysis")
    plt.tight_layout()
    plt.show()
