package com.nicksard.appsleep;

import android.content.Context;
import android.content.res.TypedArray;
import android.graphics.Typeface;
import android.graphics.drawable.GradientDrawable;
import android.util.TypedValue;
import android.view.View;
import android.view.ViewGroup;
import android.widget.LinearLayout;
import android.widget.TextView;

/** Small helpers for building the screens in code. */
final class Ui {
    private Ui() {}

    static int dp(Context c, float v) {
        return Math.round(TypedValue.applyDimension(
                TypedValue.COMPLEX_UNIT_DIP, v, c.getResources().getDisplayMetrics()));
    }

    static int themeColor(Context c, int attr) {
        TypedArray a = c.obtainStyledAttributes(new int[]{attr});
        try {
            return a.getColor(0, 0xFF888888);
        } finally {
            a.recycle();
        }
    }

    static LinearLayout column(Context c) {
        LinearLayout l = new LinearLayout(c);
        l.setOrientation(LinearLayout.VERTICAL);
        return l;
    }

    static LinearLayout row(Context c) {
        LinearLayout l = new LinearLayout(c);
        l.setOrientation(LinearLayout.HORIZONTAL);
        l.setGravity(android.view.Gravity.CENTER_VERTICAL);
        return l;
    }

    static TextView text(Context c, CharSequence s, float sp) {
        TextView t = new TextView(c);
        t.setText(s);
        t.setTextSize(TypedValue.COMPLEX_UNIT_SP, sp);
        return t;
    }

    static TextView heading(Context c, CharSequence s) {
        TextView t = text(c, s, 18);
        t.setTypeface(Typeface.DEFAULT_BOLD);
        t.setPadding(0, dp(c, 20), 0, dp(c, 8));
        return t;
    }

    static TextView secondary(Context c, CharSequence s) {
        TextView t = text(c, s, 14);
        t.setAlpha(0.75f);
        return t;
    }

    /** A rounded, subtly tinted card. */
    static LinearLayout card(Context c) {
        LinearLayout l = column(c);
        int pad = dp(c, 16);
        l.setPadding(pad, pad, pad, pad);
        GradientDrawable bg = new GradientDrawable();
        bg.setCornerRadius(dp(c, 16));
        int accent = themeColor(c, android.R.attr.colorAccent);
        bg.setColor((accent & 0x00FFFFFF) | 0x1F000000);
        l.setBackground(bg);
        LinearLayout.LayoutParams lp = new LinearLayout.LayoutParams(
                ViewGroup.LayoutParams.MATCH_PARENT, ViewGroup.LayoutParams.WRAP_CONTENT);
        lp.bottomMargin = dp(c, 12);
        l.setLayoutParams(lp);
        return l;
    }

    static LinearLayout.LayoutParams weight1() {
        return new LinearLayout.LayoutParams(0, ViewGroup.LayoutParams.WRAP_CONTENT, 1f);
    }

    static void gone(View v, boolean gone) {
        v.setVisibility(gone ? View.GONE : View.VISIBLE);
    }
}
