/* eslint-disable require-jsdoc*/
(function libraryWrapper(window) {
  function defineLibrary() {
    const VemaLib = {};

    class Vec3 {

      constructor(x, y, z) {
        this.f = [x, y, z];
      }

      length() {
        return this.calcLength();
      }

      lengthSquared() {
        const f = this.f;
        return (f[0] * f[0]) + (f[1] * f[1]) + (f[2] * f[2]);
      }

      calcLength() {
        const f = this.f;
        return Math.sqrt((f[0] * f[0]) + (f[1] * f[1]) + (f[2] * f[2]));
      }

      normalized() {
        const l = this.length;
        if (l > 0) {
          return new Vec3(this.f[0] / l, this.f[1] / l, this.f[2] / l);
        }
        return new Vec3(0, 0, 0);
      }

      addToThis(v) {
        this.f[0] += v.f[0];
        this.f[1] += v.f[1];
        this.f[2] += v.f[2];
      }

      divide(a) {
        return new Vec3(this.f[0] / a, this.f[1] / a, this.f[2] / a);
      }

      subtract(v) {
        return new Vec3(this.f[0] - v.f[0], this.f[1] - v.f[1], this.f[2] - v.f[2]);
      }

      add(v) {
        return new Vec3(this.f[0] + v.f[0], this.f[1] + v.f[1], this.f[2] + v.f[2]);
      }

      multiple(a) {
        return new Vec3(this.f[0] * a, this.f[1] * a, this.f[2] * a);
      }

      negative() {
        return new Vec3(-this.f[0], -this.f[1], -this.f[2]);
      }

      cross(v) {
        const f = this.f;
        return new Vec3((f[1] * v.f[2]) - (f[2] * v.f[1]),
                        (f[2] * v.f[0]) - (f[0] * v.f[2]),
                        (f[0] * v.f[1]) - (f[1] * v.f[0]));
      }

      dot(v) {
        const f = this.f;
        return (f[0] * v.f[0]) + (f[1] * v.f[1]) + (f[2] * v.f[2]);
      }

      clone() {
        return new Vec3(this.f[0], this.f[1], this.f[2]);
      }

      setZero() {
        this.f[0] = 0;
        this.f[1] = 0;
        this.f[2] = 0;
        return this;
      }
    }

    class Matrix4 {
      constructor() {
        this.m = [1, 0, 0, 0,
          0, 1, 0, 0,
          0, 0, 1, 0,
          0, 0, 0, 1,
        ];
      }

      makeFrustum(left, right,
        bottom, top,
        znear, zfar) {
        const X = (2 * znear) / (right - left);
        const Y = (2 * znear) / (top - bottom);
        const A = (right + left) / (right - left);
        const B = (top + bottom) / (top - bottom);
        const C = -(zfar + znear) / (zfar - znear);
        const D = (-2 * zfar * znear) / (zfar - znear);

        this.m = [X, 0, 0, 0,
          0, Y, 0, 0,
          A, B, C, -1,
          0, 0, D, 0,
        ];
        return this.m;
      }

      makePerspective(fovy, aspect, znear, zfar) {
        const ymax = znear * Math.tan((fovy * Math.PI) / 360.0);
        const ymin = -ymax;
        const xmin = ymin * aspect;
        const xmax = ymax * aspect;

        return this.makeFrustum(xmin, xmax, ymin, ymax, znear, zfar);
      }

      flatten() {
        return this.m;
      }

      setMatrixData(m) {  // toimiiko tämä???
        if (m.length === 16) {
          this.m = [
            m[0], m[1], m[2], m[3],
            m[4], m[5], m[6], m[7],
            m[8], m[9], m[10], m[11],
            m[12], m[13], m[14], m[15],
          ];
        }
      }

      multiple(m) {
        const m0 = new Matrix4();
        const m1 = m.m;
        const tm = this.m;
        m0.setMatrixData(
          [
            (tm[0] * m1[0]) + (tm[1] * m1[4]) + (tm[2] * m1[8]) + (tm[3] * m1[12]),
            (tm[0] * m1[1]) + (tm[1] * m1[5]) + (tm[2] * m1[9]) + (tm[3] * m1[13]),
            (tm[0] * m1[2]) + (tm[1] * m1[6]) + (tm[2] * m1[10]) + (tm[3] * m1[14]),
            (tm[0] * m1[3]) + (tm[1] * m1[7]) + (tm[2] * m1[11]) + (tm[3] * m1[15]),

            (tm[4] * m1[0]) + (tm[5] * m1[4]) + (tm[6] * m1[8]) + (tm[7] * m1[12]),
            (tm[4] * m1[1]) + (tm[5] * m1[5]) + (tm[6] * m1[9]) + (tm[7] * m1[13]),
            (tm[4] * m1[2]) + (tm[5] * m1[6]) + (tm[6] * m1[10]) + (tm[7] * m1[14]),
            (tm[4] * m1[3]) + (tm[5] * m1[7]) + (tm[6] * m1[11]) + (tm[7] * m1[15]),

            (tm[8] * m1[0]) + (tm[9] * m1[4]) + (tm[10] * m1[8]) + (tm[11] * m1[12]),
            (tm[8] * m1[1]) + (tm[9] * m1[5]) + (tm[10] * m1[9]) + (tm[11] * m1[13]),
            (tm[8] * m1[2]) + (tm[9] * m1[6]) + (tm[10] * m1[10]) + (tm[11] * m1[14]),
            (tm[8] * m1[3]) + (tm[9] * m1[7]) + (tm[10] * m1[11]) + (tm[11] * m1[15]),

            (tm[12] * m1[0]) + (tm[13] * m1[4]) + (tm[14] * m1[8]) + (tm[15] * m1[12]),
            (tm[12] * m1[1]) + (tm[13] * m1[5]) + (tm[14] * m1[9]) + (tm[15] * m1[13]),
            (tm[12] * m1[2]) + (tm[13] * m1[6]) + (tm[14] * m1[10]) + (tm[15] * m1[14]),
            (tm[12] * m1[3]) + (tm[13] * m1[7]) + (tm[14] * m1[11]) + (tm[15] * m1[15]),
          ]
        );

        return m0;
      }


      multipleV(v) {
        const m = this.m;
        return new Vec3((v.f[0] * m[0]) + (v.f[1] * m[4]) + (v.f[2] * m[8]) + m[12],
          (v.f[0] * m[1]) + (v.f[1] * m[5]) + (v.f[2] * m[9]) + m[13],
          (v.f[0] * m[2]) + (v.f[1] * m[6]) + (v.f[2] * m[10]) + m[14]);
      }

      rotationYMatrix(t) {
        const c = Math.cos(t);
        const s = Math.sin(t);
        this.m = [c, 0, -s, 0,
          0, 1, 0, 0,
          s, 0, c, 0,
          0, 0, 0, 1,
        ];
        return this.m;
      }

      rotationXMatrix(t) {
        const c = Math.cos(t);
        const s = Math.sin(t);
        this.m = [1, 0, 0, 0,
          0, c, -s, 0,
          0, s, c, 0,
          0, 0, 0, 1,
        ];
        return this.m;
      }

      rotationZMatrix(t) {
        const c = Math.cos(t);
        const s = Math.sin(t);
        this.m = [c, -s, 0, 0,
          s, c, 0, 0,
          0, 0, 1, 0,
          0, 0, 0, 1,
        ];
        return this.m;
      }

      rotationQMatrix(q) {
        this.m = [1 - (2 * q[1] * q[1]) - (2 * q[2] * q[2]),
          (2 * q[0] * q[1]) + (2 * q[2] * q[3]),
          (2 * q[0] * q[2]) - (2 * q[1] * q[3]),
          0,
          (2 * q[0] * q[1]) - (2 * q[2] * q[3]),
          1 - (2 * q[0] * q[0]) - (2 * q[2] * q[2]),
          (2 * q[1] * q[2]) + (2 * q[0] * q[3]),
          0,
          (2 * q[0] * q[2]) + (2 * q[1] * q[3]),
          (2 * q[1] * q[2]) - (2 * q[0] * q[3]),
          1 - (2 * q[0] * q[0]) - (2 * q[1] * q[1]),
          0,
          0, 0, 0, 1,
        ];
        return this.m;
      }

      rotateMatrix(theta, normal) {
        const c = Math.cos(theta);
        const s = Math.sin(theta);
        const bigC = 1 - c;
        const f = normal.f;

        this.m = [
          (f[0] * f[0] * bigC) + c,
          (f[0] * f[1] * bigC) - (f[2] * s),
          (f[0] * f[2] * bigC) + (f[1] * s),
          0,
          (f[1] * f[0] * bigC) + (f[2] * s),
          (f[1] * f[1] * bigC) + c,
          (f[1] * f[2] * bigC) - (f[0] * s),
          0,
          (f[2] * f[0] * bigC) - (f[1] * s),
          (f[2] * f[1] * bigC) + (f[0] * s),
          (f[2] * f[2] * bigC) + c,
          0,
          0, 0, 0, 1,
        ];
        return this.m;
      }

      translate(x, y, z) {
        this.m = [
          1, 0, 0, 0,
          0, 1, 0, 0,
          0, 0, 1, 0,
          x, y, z, 1,
        ];
        return this.m;
      }

      scale(x, y, z) {
        this.m = [
          x, 0, 0, 0,
          0, y, 0, 0,
          0, 0, z, 0,
          0, 0, 0, 1
        ];
        return this.m;
      }

      randomRotationMatrix() {
        const theta = Math.random() * Math.PI * 2;
        const phi = Math.random() * Math.PI * 2;
        const z = Math.random() * 2.0;

        const r = Math.sqrt(z);
        const Vx = Math.sin(phi) * r;
        const Vy = Math.cos(phi) * r;
        const Vz = Math.sqrt(2.0 - z);

        const st = Math.sin(theta);
        const ct = Math.cos(theta);
        const Sx = (Vx * ct) - (Vy * st);
        const Sy = (Vx * st) + (Vy * ct);

        this.m = [
          (Vx * Sx) - ct, (Vx * Sy) - st, Vx * Vz, 0,
          (Vy * Sx) + st, (Vy * Sy) - ct, Vy * Vz, 0,
          Vz * Sx, Vz * Sy, 1.0 - z, 0,
          0, 0, 0, 1,
        ];
        return this.m;
      }

      makeLookAt(direction, up, right, position) {
        const p = new Matrix4();
        p.setMatrixData(
          [1, 0, 0, 0,
            0, 1, 0, 0,
            0, 0, 1, 0,
            -position.f[0], -position.f[1], -position.f[2], 1]);

        const r = new Matrix4();
        r.setMatrixData(
          [right.f[0], up.f[0], direction.f[0], 0,
            right.f[1], up.f[1], direction.f[1], 0,
            right.f[2], up.f[2], direction.f[2], 0,
            0, 0, 0, 1]);
        const result = p.multiple(r);
        this.m = result.m;
        return this.m;
      }

      makeLookAt2(direction, up, right, position) {
        const p = new Matrix4();
        p.setMatrixData(
          [1, 0, 0, 0,
            0, 1, 0, 0,
            0, 0, 1, 0,
            -position.f[0], -position.f[1], -position.f[2], 1]);

        const r = new Matrix4();
        r.setMatrixData(
          [right.f[0], right.f[1], right.f[2], 0,
            up.f[0], up.f[1], up.f[2], 0,
            direction.f[0], direction.f[1], direction.f[2], 0,
            0, 0, 0, 1]);
        const result = p.multiple(r);
        this.m = result.m;
        return this.m;
      }


      invert() {
        const inv = [1, 0, 0, 0,
          0, 1, 0, 0,
          0, 0, 1, 0,
          0, 0, 0, 1];
        let det;
        const m1 = this.m;

        inv[0] = m1[5] * m1[10] * m1[15] -
             m1[5] * m1[11] * m1[14] -
             m1[9] * m1[6] * m1[15] +
             m1[9] * m1[7] * m1[14] +
             m1[13] * m1[6] * m1[11] -
             m1[13] * m1[7] * m1[10];

        inv[4] = -m1[4] * m1[10] * m1[15] +
              m1[4] * m1[11] * m1[14] +
              m1[8] * m1[6] * m1[15] -
              m1[8] * m1[7] * m1[14] -
              m1[12] * m1[6] * m1[11] +
              m1[12] * m1[7] * m1[10];

        inv[8] = m1[4] * m1[9] * m1[15] -
             m1[4] * m1[11] * m1[13] -
             m1[8] * m1[5] * m1[15] +
             m1[8] * m1[7] * m1[13] +
             m1[12] * m1[5] * m1[11] -
             m1[12] * m1[7] * m1[9];

        inv[12] = -m1[4] * m1[9] * m1[14] +
               m1[4] * m1[10] * m1[13] +
               m1[8] * m1[5] * m1[14] -
               m1[8] * m1[6] * m1[13] -
               m1[12] * m1[5] * m1[10] +
               m1[12] * m1[6] * m1[9];

        inv[1] = -m1[1] * m1[10] * m1[15] +
              m1[1] * m1[11] * m1[14] +
              m1[9] * m1[2] * m1[15] -
              m1[9] * m1[3] * m1[14] -
              m1[13] * m1[2] * m1[11] +
              m1[13] * m1[3] * m1[10];

        inv[5] = m1[0] * m1[10] * m1[15] -
             m1[0] * m1[11] * m1[14] -
             m1[8] * m1[2] * m1[15] +
             m1[8] * m1[3] * m1[14] +
             m1[12] * m1[2] * m1[11] -
             m1[12] * m1[3] * m1[10];

        inv[9] = -m1[0] * m1[9] * m1[15] +
              m1[0] * m1[11] * m1[13] +
              m1[8] * m1[1] * m1[15] -
              m1[8] * m1[3] * m1[13] -
              m1[12] * m1[1] * m1[11] +
              m1[12] * m1[3] * m1[9];

        inv[13] = m1[0] * m1[9] * m1[14] -
              m1[0] * m1[10] * m1[13] -
              m1[8] * m1[1] * m1[14] +
              m1[8] * m1[2] * m1[13] +
              m1[12] * m1[1] * m1[10] -
              m1[12] * m1[2] * m1[9];

        inv[2] = m1[1] * m1[6] * m1[15] -
             m1[1] * m1[7] * m1[14] -
             m1[5] * m1[2] * m1[15] +
             m1[5] * m1[3] * m1[14] +
             m1[13] * m1[2] * m1[7] -
             m1[13] * m1[3] * m1[6];

        inv[6] = -m1[0] * m1[6] * m1[15] +
              m1[0] * m1[7] * m1[14] +
              m1[4] * m1[2] * m1[15] -
              m1[4] * m1[3] * m1[14] -
              m1[12] * m1[2] * m1[7] +
              m1[12] * m1[3] * m1[6];

        inv[10] = m1[0] * m1[5] * m1[15] -
              m1[0] * m1[7] * m1[13] -
              m1[4] * m1[1] * m1[15] +
              m1[4] * m1[3] * m1[13] +
              m1[12] * m1[1] * m1[7] -
              m1[12] * m1[3] * m1[5];

        inv[14] = -m1[0] * m1[5] * m1[14] +
               m1[0] * m1[6] * m1[13] +
               m1[4] * m1[1] * m1[14] -
               m1[4] * m1[2] * m1[13] -
               m1[12] * m1[1] * m1[6] +
               m1[12] * m1[2] * m1[5];

        inv[3] = -m1[1] * m1[6] * m1[11] +
              m1[1] * m1[7] * m1[10] +
              m1[5] * m1[2] * m1[11] -
              m1[5] * m1[3] * m1[10] -
              m1[9] * m1[2] * m1[7] +
              m1[9] * m1[3] * m1[6];

        inv[7] = m1[0] * m1[6] * m1[11] -
             m1[0] * m1[7] * m1[10] -
             m1[4] * m1[2] * m1[11] +
             m1[4] * m1[3] * m1[10] +
             m1[8] * m1[2] * m1[7] -
             m1[8] * m1[3] * m1[6];

        inv[11] = -m1[0] * m1[5] * m1[11] +
               m1[0] * m1[7] * m1[9] +
               m1[4] * m1[1] * m1[11] -
               m1[4] * m1[3] * m1[9] -
               m1[8] * m1[1] * m1[7] +
               m1[8] * m1[3] * m1[5];

        inv[15] = m1[0] * m1[5] * m1[10] -
              m1[0] * m1[6] * m1[9] -
              m1[4] * m1[1] * m1[10] +
              m1[4] * m1[2] * m1[9] +
              m1[8] * m1[1] * m1[6] -
              m1[8] * m1[2] * m1[5];

        // console.log(inv);
        det = m1[0] * inv[0] + m1[1] * inv[4] + m1[2] * inv[8] + m1[3] * inv[12];

        if (det == 0) {
          return new Matrix4();
        }

        det = 1.0 / det;


        const res = new Matrix4();
        res.setMatrixData([
          inv[0] * det, inv[1] * det, inv[2] * det, inv[3] * det,
          inv[4] * det, inv[5] * det, inv[6] * det, inv[7] * det,
          inv[8] * det, inv[9] * det, inv[10] * det, inv[11] * det,
          inv[12] * det, inv[13] * det, inv[14] * det, inv[15] * det]);
        return res;
      }

      transpose() {
        const m1 = this.m;
        const res = new Matrix4();
        res.setMatrixData([m1[0], m1[4], m1[8], m1[12],
          m1[1], m1[5], m1[9], m1[13],
          m1[2], m1[6], m1[10], m1[14],
          m1[3], m1[7], m1[11], m1[15]]);
        return res;
      }


    }

    VemaLib.Vec3 = Vec3;
    VemaLib.Matrix4 = Matrix4;

    return VemaLib;
  }
  if (typeof (VemaLib) === 'undefined') window.VemaLib = defineLibrary(); // eslint-disable-line no-param-reassign, no-undef
  else console.log('Library already defined.'); // eslint-disable-line no-console
}(window)); // eslint-disable-line no-undef
