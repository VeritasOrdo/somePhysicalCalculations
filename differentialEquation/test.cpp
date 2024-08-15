#include <stdio.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_errno.h>

// 定义微分方程组
int func(double t, const double y[], double f[], void *params) {
    f[0] = y[1];
    f[1] = -y[0];
    return GSL_SUCCESS;
}

// 定义雅可比矩阵（对于简单的例子可以不使用）
int jac(double t, const double y[], double *dfdy, double dfdt[], void *params) {
    gsl_matrix_view dfdy_mat = gsl_matrix_view_array(dfdy, 2, 2);
    gsl_matrix *m = &dfdy_mat.matrix;
    gsl_matrix_set(m, 0, 0, 0.0);
    gsl_matrix_set(m, 0, 1, 1.0);
    gsl_matrix_set(m, 1, 0, -1.0);
    gsl_matrix_set(m, 1, 1, 0.0);
    dfdt[0] = 0.0;
    dfdt[1] = 0.0;
    return GSL_SUCCESS;
}

int main() {
    // 定义 ODE 系统
    gsl_odeiv2_system sys = {func, jac, 2, NULL};

    // 创建步进器和控制器
    gsl_odeiv2_driver *d = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rk8pd, 1e-6, 1e-6, 0.0);

    // 设置初始条件
    double t = 0.0;
    double y[2] = {1.0, 0.0}; // y(0) = 1, y'(0) = 0

    // 进行积分
    for (int i = 1; i <= 100; i++) {
        double ti = i * 0.1;
        int status = gsl_odeiv2_driver_apply(d, &t, ti, y);
        if (status != GSL_SUCCESS) {
            printf("error: %d\n", status);
            break;
        }
        printf("t = %.2f, y[0] = %.5f, y[1] = %.5f\n", t, y[0], y[1]);
    }

    // 释放资源
    gsl_odeiv2_driver_free(d);
    return 0;
}