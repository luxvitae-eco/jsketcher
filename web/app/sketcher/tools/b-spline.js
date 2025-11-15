import { Tool } from "./tool";
import { Segment } from "../shapes/segment";
import { EndPoint } from "../shapes/point";
import { IBSplineOpts, CentripetalParameterMethod, CPointsCalculator, BSpline } from "../shapes/b-spline";
import Vector from "math/vector";
import { TOLERANCE, arePointsEqual } from "math/equality";

export class BSplineTool extends Tool {
  constructor(viewer) {
    super("basic spline curve", viewer);
    this.init();
    this._v = new Vector();
  }

  init() {
    this.degree = 3;
    this.fPoints = [];
    this.curve = null;
    this.otherCurveEndPoint = null;
  }

  restart() {
    this.init();
    this.sendHint("specify first point");
  }

  cleanup(e) {
    this.viewer.cleanSnap();
  }

  mouseup(e) {
    const p = this.viewer.screenToModel(e);
    const length = this.fPoints.length;
    if (length && arePointsEqual(this.fPoints[length - 1], p)) {
      return;
    }
    const point = new EndPoint(p.x, p.y);
    this.fPoints.push(point);
    if (this.fPoints.length < 2) {
      this.curve = new Segment(this.fPoints[0].x, this.fPoints[0].y, point.x, point.y);
      this.viewer.activeLayer.add(this.curve);
    } else if (this.fPoints.length == 2) {
      const opts = {
        degree: this.degree,
        cPoints: [this.fPoints[0], this.fPoints[0], this.fPoints[1], this.fPoints[1], this.fPoints[1]],
        fPoints: [...this.fPoints, new EndPoint(p.x + 0.1, p.y + 0.1)],
        kValues: [0, 0, 0, 1, 1, 1, 0, 0, 0],
      };
      this.viewer.activeLayer.remove(this.curve);
      this.curve = new BSpline(opts, true, false);
      this.curve.update();
      this.viewer.activeLayer.add(this.curve);
    } else {
      this.curve.removeFPoint();
      this.curve.addFPoint(point);
      this.curve.update();
    }
    if (this.curve !== null) {
      this.curve.stabilize(this.viewer);
    }
    this.viewer.refresh();
  }

  mousemove(e) {
    if (this.curve == null) {
      return;
    }
    const p = this.viewer.screenToModel(e);
    if (this.fPoints.length < 2) {
      this.curve.b.x = p.x;
      this.curve.b.y = p.y;
    } else {
      if (this.curve.fPoints.length == this.fPoints.length) {
        const point = new EndPoint(p.x, p.y);
        this.curve.addFPoint(point);
      } else {
        this.curve.fPoints[this.fPoints.length].x = p.x;
        this.curve.fPoints[this.fPoints.length].y = p.y;
      }
      this.curve.update();
    }

    this.viewer.refresh();
  }
}
