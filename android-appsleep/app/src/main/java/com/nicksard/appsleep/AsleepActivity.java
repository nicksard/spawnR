package com.nicksard.appsleep;

import android.app.Activity;
import android.content.Intent;
import android.content.pm.ApplicationInfo;
import android.content.pm.PackageManager;
import android.graphics.drawable.Drawable;
import android.os.Bundle;
import android.view.Gravity;
import android.view.ViewGroup;
import android.widget.Button;
import android.widget.ImageView;
import android.widget.LinearLayout;
import android.widget.TextView;

/** The screen shown after a sleeping app was closed. */
public class AsleepActivity extends Activity {
    static final String EXTRA_PACKAGE = "pkg";
    static final String EXTRA_SCHEDULE_ID = "schedule";

    @Override
    protected void onCreate(Bundle savedInstanceState) {
        super.onCreate(savedInstanceState);
        String pkg = getIntent().getStringExtra(EXTRA_PACKAGE);
        Schedule s = ScheduleStore.find(this, getIntent().getStringExtra(EXTRA_SCHEDULE_ID));

        PackageManager pm = getPackageManager();
        CharSequence label = pkg;
        Drawable icon = null;
        try {
            ApplicationInfo ai = pm.getApplicationInfo(pkg, 0);
            label = pm.getApplicationLabel(ai);
            icon = pm.getApplicationIcon(ai);
        } catch (PackageManager.NameNotFoundException | NullPointerException ignored) {
        }

        LinearLayout root = Ui.column(this);
        root.setGravity(Gravity.CENTER);
        int pad = Ui.dp(this, 32);
        root.setPadding(pad, pad, pad, pad);

        if (icon != null) {
            ImageView iv = new ImageView(this);
            iv.setImageDrawable(icon);
            iv.setAlpha(0.5f);
            int size = Ui.dp(this, 72);
            root.addView(iv, new LinearLayout.LayoutParams(size, size));
        }

        TextView moon = Ui.text(this, "🌙", 40);
        moon.setGravity(Gravity.CENTER);
        moon.setPadding(0, Ui.dp(this, 12), 0, 0);
        root.addView(moon);

        TextView title = Ui.text(this, label + " is asleep", 24);
        title.setGravity(Gravity.CENTER);
        title.setPadding(0, Ui.dp(this, 8), 0, Ui.dp(this, 8));
        root.addView(title);

        String detail = s == null
                ? "It's in a sleep window right now."
                : "It wakes up at " + Schedule.formatMinute(s.endMinute)
                    + ".\nSchedule: " + s.name;
        TextView body = Ui.secondary(this, detail);
        body.setGravity(Gravity.CENTER);
        root.addView(body);

        Button ok = new Button(this);
        ok.setText("OK");
        ok.setOnClickListener(v -> goHome());
        LinearLayout.LayoutParams lp = new LinearLayout.LayoutParams(
                ViewGroup.LayoutParams.WRAP_CONTENT, ViewGroup.LayoutParams.WRAP_CONTENT);
        lp.topMargin = Ui.dp(this, 24);
        root.addView(ok, lp);

        setContentView(root);
    }

    @Override
    public void onBackPressed() {
        goHome();
    }

    private void goHome() {
        startActivity(new Intent(Intent.ACTION_MAIN)
                .addCategory(Intent.CATEGORY_HOME)
                .addFlags(Intent.FLAG_ACTIVITY_NEW_TASK));
        finish();
    }
}
